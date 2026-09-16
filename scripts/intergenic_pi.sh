download 
Orthogroups.tsv
SCO.txt
orthodiver_agg

What do i need to run dfe alpha?
Find other notebook to make plots

Definitely make sure NAs and 0s are correct in gene pop

# construct divergence_file.txt from *{0,4}d_sites_summary.txt 
#1,total number of 0d sites, number of 0d fixed differences
#2,number of 4d sites, number of 4d fixed differences

#afusca v dicmin
1 6588078 253571 
0 1345742 813059

#afusca v liplubb
1 5339072 181148
0 1130462 599042

#afusca v smivir
1 10187118 316453
0 2163103 965776

data/workdir/ortholog_pop_gen/75/2e04606d0dcbcf6a2dbf46e04e73e7/WW5-5.ref.results.0d_sites_summary.txt
data/workdir/ortholog_pop_gen/75/2e04606d0dcbcf6a2dbf46e04e73e7/WW5-5.ref.results.4d_sites_summary.txt

SVIR A
lambda 0.673940 selected_divergence 0.030700 alpha -7.277078 omega_A -0.331497

SVIR X
lambda 0.675593 selected_divergence 0.030921 alpha -3.763470 omega_A -0.172252

main things:
resampling genes for dfe alpha sfs
alpha from dnds analysis
chromosome comparisons
sex-biased gene slices
make sure to get paper extension and reply to emails i need to

Neut A dfe
N1 100 N2 4 t2 6.3716 N3 50 t3 15.7455 Nw 47.68 f0 0.994360531 L -144945.9450

Neut X dfe
N1 100 N2 13 t2 92.7581 Nw 15.46 f0 0.997730359 L -17715.4949

Neut A egf
0 8
0 0.9936779804
1 0.0014183859
2 0.0005515046
3 0.0003647233
4 0.0002541172
5 0.0001904317
6 0.0001500978
7 0.0001224329
8 0.0001025353

Weird?
Neut X egf
0 26
0 0.9976265381
1 0.0005240675
2 0.0002301599
3 0.0001652380
4 0.0001260489
5 0.0001028730
6 0.0000876315
7 0.0000767376
8 0.0000685514
9 0.0000621802
10 0.0000570821
11 0.0000529094
12 0.0000494302
13 0.0000464841
14 0.0000439567
15 0.0000417639
16 0.0000398427
17 0.0000381450
18 0.0000366332
19 0.0000352776
20 0.0000340544
21 0.0000329442
22 0.0000319311
23 0.0000310019
24 0.0000301452
25 0.0000293515
26 0.0000286122

Sel A dfe
N1 100 N2 4 t2 6.3716 N3 50 t3 15.7455 Nw 47.68 b 0.1305 Es -4.286763 f0 0.994360531 L -99172.8867

Sel X dfe
N1 100 N2 4 t2 6.3716 N3 50 t3 15.7455 Nw 47.68 b 0.2179 Es -4.373423 f0 0.994360531 L -15317.1044

Sel A egf
1 8
0 0.9970165750
1 0.0009633369
2 0.0003083849
3 0.0001914899
4 0.0001270669
5 0.0000917327
6 0.0000701990
7 0.0000559149
8 0.0000459505

Sel X egf
1 8
0 0.9977682809
1 0.0008592400
2 0.0002453520
3 0.0001470115
4 0.0000949126
5 0.0000670698
6 0.0000504232
7 0.0000395460
8 0.0000320506

prop mut in s ranges

A
n2 47.680000
es -4.286763
beta 0.130500
mean_ns 204.392860
lower 0.000000 upper 1.000000 area 0.407405
lower 1.000000 upper 10.000000 area 0.142438
lower 10.000000 upper 100.000000 area 0.187890
lower 100.000000 upper -99.000000 area 0.262267
Total area 1.000000


X
n2 47.680000
es -4.373423
beta 0.217900
mean_ns 208.524809
lower 0.000000 upper 1.000000 area 0.245255
lower 1.000000 upper 10.000000 area 0.159124
lower 10.000000 upper 100.000000 area 0.252574
lower 100.000000 upper -99.000000 area 0.343047
Total area 1.000000

X 0D
0.00140 [0.00124,0.00156]
X 4D
0.00336 [0.00285,0.00388]
A 0D
0.00417 [0.00406,0.00427]
A 4D
0.00491 [0.00474,0.00508]


parallel -j1 'grep {} allacma_fusca.callable.intergenic.rm_repeats.bed > {}_intergenic.bed' :::: chroms.txt
parallel './calculatePiBed.py -v allacma_fusca.hard_filtered.sorted.vcf.gz -b {}_intergenic.bed -g allacma_fusca.genomefile -n {} -l intergenic -o {}.intergenic.pi.tsv' :::: chroms.txt
cat *1.intergenic.pi.tsv > tmp
grep -v intergenic tmp > all_chrs.intergenic.pi.tsv

check for pi vs window size correlation
easySFS using autosome and X slices

check dont need to subtract repeats from the callable intergenic bed
is vcftools counting over the whole window or
