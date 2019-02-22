#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import system
from meld import comm, vault
from meld import parse
from meld.system import patchers


N_REPLICAS = 2
N_STEPS = 10000
BLOCK_SIZE = 100


def gen_state(s, index):
    pos = s._coordinates
    pos = pos - np.mean(pos, axis=0)
    box = [999., 999., 999.]
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy, box)


def setup_system():
    # load the sequence
    # sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    # n_res = len(sequence.split())

    # setup the patchers
    patcher = patchers.VirtualSpinLabel(
        {
            122: 'OND',
            37: 'OND',
            48: 'OND',
            138: 'OND',
            81: 'OND',
            12: 'OND',
            105: 'OND',
            112: 'OND',
            29: 'OND'})


    # build the system
    p = system.ProteinMoleculeFromPdbFile('cam2smmlck.pdb')
    b = system.SystemBuilder()
    s = b.build_system_from_molecules([p], [patcher])
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.5, 300., 550.)

    # add restraint
    # scaler = s.restraints.create_scaler('constant')
    # r = s.restraints.create_restraint('distance',
    #                                   scaler,
    #                                   r1=0.0,
    #                                   r2=0.0,
    #                                   r3=.45,
    #                                   r4=.65,
    #                                   k=250.,
    #                                   atom_1_res_index=1,
    #                                   atom_1_name='OND',
    #                                   atom_2_res_index=16,
    #                                   atom_2_name='OND')
    # s.restraints.add_as_always_active(r)

    # create the options
    options = system.RunOptions(solvation='implicit')
    options.implicit_solvent_model = 'gbNeck'
    options.use_big_timestep = False
    options.cutoff = 1.8

    options.use_amap = True
    options.amap_beta_bias = 1.0
    options.timesteps = 1428
    options.minimize_steps = 100

    # create a store
    store = vault.DataStore(s.n_atoms,
                            N_REPLICAS,
                            s.get_pdb_writer(),
                            block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS,
                                       adaptation_policy=policy_1)

    remd_runner = master_runner.MasterReplicaExchangeRunner(
        N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms

setup_system()
