@time Qin_t, Qout_t, A_t = modes_init_conds(Ne-2, Ne-1, 0, fps, Tf)
@time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

