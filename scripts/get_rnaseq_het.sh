VCF=data/results/var_call/allacma_fusca_rnaseq/allacma_fusca.hard_filtered.sorted.vcf.gz

bcftools view \
    -r OX359250.1,OX359249.1 \
    "$VCF" \
| bcftools query -f '[%SAMPLE\t%GT\n]' \
| awk '
{
    gt=$2
    if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0")
        het[$1]++
}
END{
    print "sample\tsex\ttype\thet_sites"
    for(s in het){
        split(s,a,"_")
        sex=a[2]
        n=a[3]+0
        type=(n<=5 ? "pooled" : "single")
        print s"\t"sex"\t"type"\t"het[s]
    }
}' | sort -k2,2 -k3,3 -k1,1
