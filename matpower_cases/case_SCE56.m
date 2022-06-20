function mpc = case_SCE56
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	1	1	100	0;
	2	1	0	0	0	0	1	1	0	1	1	100	0;
	3	1	0	0	0	0	1	1	0	1	1	100	0;
    4	1	0	0	0	0	1	1	0	1	1	100	0;
    5	1	0	0	0	0	1	1	0	1	1	100	0;
    6	1	0	0	0	0	1	1	0	1	1	100	0;
    7	1	0	0	0	0	1	1	0	1	1	100	0;
    8	1	0	0	0	0	1	1	0	1	1	100	0;
    9	1	0	0	0	0	1	1	0	1	1	100	0;
    10	1	0	0	0	0	1	1	0	1	1	100	0;
    11	1	0	0	0	0	1	1	0	1	1	100	0;
    12	1	0	0	0	0	1	1	0	1	1	100	0;
    13	1	0	0	0	0	1	1	0	1	1	100	0;
    14	1	0	0	0	0	1	1	0	1	1	100	0;
    15	1	0	0	0	0	1	1	0	1	1	100	0;
    16	1	0	0	0	0	1	1	0	1	1	100	0;
    17	1	0	0	0	0	1	1	0	1	1	100	0;
    18	1	0	0	0	0	1	1	0	1	1	100	0;
    19	1	0	0	0	0	1	1	0	1	1	100	0;
    20	1	0	0	0	0	1	1	0	1	1	100	0;
    21	1	0	0	0	0	1	1	0	1	1	100	0;
    22	1	0	0	0	0	1	1	0	1	1	100	0;
    23	1	0	0	0	0	1	1	0	1	1	100	0;
    24	1	0	0	0	0	1	1	0	1	1	100	0;
    25	1	0	0	0	0	1	1	0	1	1	100	0;
    26	1	0	0	0	0	1	1	0	1	1	100	0;
    27	1	0	0	0	0	1	1	0	1	1	100	0;
    28	1	0	0	0	0	1	1	0	1	1	100	0;
    29	1	0	0	0	0	1	1	0	1	1	100	0;
    30	1	0	0	0	0	1	1	0	1	1	100	0;
    31	1	0	0	0	0	1	1	0	1	1	100	0;
    32	1	0	0	0	0	1	1	0	1	1	100	0;
    33	1	0	0	0	0	1	1	0	1	1	100	0;
    34	1	0	0	0	0	1	1	0	1	1	100	0;
    35	1	0	0	0	0	1	1	0	1	1	100	0;
    36	1	0	0	0	0	1	1	0	1	1	100	0;
    37	1	0	0	0	0	1	1	0	1	1	100	0;
    38	1	0	0	0	0	1	1	0	1	1	100	0;
    39	1	0	0	0	0	1	1	0	1	1	100	0;
    40	1	0	0	0	0	1	1	0	1	1	100	0;
    41	1	0	0	0	0	1	1	0	1	1	100	0;
    42	1	0	0	0	0	1	1	0	1	1	100	0;
    43	1	0	0	0	0	1	1	0	1	1	100	0;
    44	1	0	0	0	0	1	1	0	1	1	100	0;
    45	1	0	0	0	0	1	1	0	1	1	100	0;
    46	1	0	0	0	0	1	1	0	1	1	100	0;
    47	1	0	0	0	0	1	1	0	1	1	100	0;
    48	1	0	0	0	0	1	1	0	1	1	100	0;
    49	1	0	0	0	0	1	1	0	1	1	100	0;
    50	1	0	0	0	0	1	1	0	1	1	100	0;
    51	1	0	0	0	0	1	1	0	1	1	100	0;
    52	1	0	0	0	0	1	1	0	1	1	100	0;
    53	1	0	0	0	0	1	1	0	1	1	100	0;
    54	1	0	0	0	0	1	1	0	1	1	100	0;
    55	1	0	0	0	0	1	1	0	1	1	100	0;
    56	1	0	0	0	0	1	1	0	1	1	100	0;
    
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1	0	0	100	-100	1	100	1	100	-100	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1	2	0.160	0.388	0	0	0	0	0	0	1	0	0;
	2	3	0.824	0.315	0	0	0	0	0	0	1	0	0;
    2	4	0.144	0.349	0	0	0	0	0	0	1	0	0;
    4	5	1.026	0.421	0	0	0	0	0	0	1	0	0;
    4	6	0.741	0.466	0	0	0	0	0	0	1	0	0;
    4	7	0.528	0.468	0	0	0	0	0	0	1	0	0;
    4	20	0.138	0.334	0	0	0	0	0	0	1	0	0;
    7	8	0.358	0.314	0	0	0	0	0	0	1	0	0;
    8	9	2.032	0.798	0	0	0	0	0	0	1	0	0;
    8	10	0.502	0.441	0	0	0	0	0	0	1	0	0;
    10	11	0.372	0.327	0	0	0	0	0	0	1	0	0;
    11	12	1.431	0.999	0	0	0	0	0	0	1	0	0;
    11	13	0.429	0.377	0	0	0	0	0	0	1	0	0;
    13	14	0.671	0.257	0	0	0	0	0	0	1	0	0;
    13	15	0.457	0.401	0	0	0	0	0	0	1	0	0;
    15	16	1.008	0.385	0	0	0	0	0	0	1	0	0;
    15	17	0.153	0.134	0	0	0	0	0	0	1	0	0;
    17	18	0.971	0.722	0	0	0	0	0	0	1	0	0;
    18	19	1.885   0.721	0	0	0	0	0	0	1	0	0;
    20	21	1.818   0.695	0	0	0	0	0	0	1	0	0;
    20	23	0.225   0.542	0	0	0	0	0	0	1	0	0;
    21	22	1.818   0.695	0	0	0	0	0	0	1	0	0;
    23	24	0.127   0.028	0	0	0	0	0	0	1	0	0;
    23	25	0.284   0.687	0	0	0	0	0	0	1	0	0;
    25	26	0.171   0.414	0	0	0	0	0	0	1	0	0;
    26	27	0.414   0.386	0	0	0	0	0	0	1	0	0;
    26	32	0.205   0.495	0	0	0	0	0	0	1	0	0;
    27	28	0.210   0.196	0	0	0	0	0	0	1	0	0;
    28	29	0.395   0.369	0	0	0	0	0	0	1	0	0;
    29	30	0.248   0.232	0	0	0	0	0	0	1	0	0;
    30	31	0.279   0.260	0	0	0	0	0	0	1	0	0;
    32	33	0.263   0.073	0	0	0	0	0	0	1	0	0;
    32	34	0.071   0.171	0	0	0	0	0	0	1	0	0;
    34	35	0.625   0.273	0	0	0	0	0	0	1	0	0;
    34	36	0.510   0.209	0	0	0	0	0	0	1	0	0;
    34	38	1.062   0.406	0	0	0	0	0	0	1	0	0;
    34	41	0.115	0.278	0	0	0	0	0	0	1	0	0;
    36	37	2.018   0.829	0	0	0	0	0	0	1	0	0;
    38	39	0.610   0.238	0	0	0	0	0	0	1	0	0;
    39	40	2.349	0.964	0	0	0	0	0	0	1	0	0;
    41	42	0.159	0.384	0	0	0	0	0	0	1	0	0;
    41	47	0.157	0.379	0	0	0	0	0	0	1	0	0;
    42	43	0.934	0.383	0	0	0	0	0	0	1	0	0;
    42	44	0.506	0.163	0	0	0	0	0	0	1	0	0;
    42	45	0.095	0.195	0	0	0	0	0	0	1	0	0;
    42	46	1.915	0.769	0	0	0	0	0	0	1	0	0;
    47	48	1.641	0.670	0	0	0	0	0	0	1	0	0;
    47	49	0.081	0.196	0	0	0	0	0	0	1	0	0;
    49	50	1.727	0.709	0	0	0	0	0	0	1	0	0;
    49	51	0.112	0.270	0	0	0	0	0	0	1	0	0;
    51	52	0.674	0.275	0	0	0	0	0	0	1	0	0;
    51	53	0.070	0.170	0	0	0	0	0	0	1	0	0;
    53	54	2.041	0.780	0	0	0	0	0	0	1	0	0;
    53	55	0.813	0.334	0	0	0	0	0	0	1	0	0;
    53	56	0.141	0.340	0	0	0	0	0	0	1	0	0;
];