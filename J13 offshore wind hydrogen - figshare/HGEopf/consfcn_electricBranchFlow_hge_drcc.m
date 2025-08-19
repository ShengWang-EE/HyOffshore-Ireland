function h = consfcn_electricBranchFlow_hge_drcc(Va, phi_Va,u_Va, v_Va, mpc, il, sigma, epsilon)
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
    upf = mpc.branch(il, RATE_A) / mpc.baseMVA - Pfinj(il);
    upt = mpc.branch(il, RATE_A) / mpc.baseMVA + Pfinj(il);
%% unpack data
h1 = Bf(il,:)*Va - upf;
h2 = -upt - Bf(il,:)*Va;
h3 = u_Va.^2 + sigma .* (Bf(il,:)*phi_Va).^2 - (upt-v_Va).^2 * epsilon;
h4 = -(u_Va+v_Va) - (Bf(il,:)*Va);
h5 = (Bf(il,:)*Va) - (u_Va+v_Va);
h6 = -u_Va;
h7 = -v_Va;
h8 = v_Va - upt;
h = [h1;h2;h3;h4;h5;h6;h7;h8];
end