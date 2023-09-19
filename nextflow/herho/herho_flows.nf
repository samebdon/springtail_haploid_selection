include {herhoTally; herhoStandAlone} from './herho_tasks.nf'

workflow herho_flow {
	take:
		vcf
		vcf_index
		bed
		species
	main:
		herhoTally(vcf, vcf_index, bed, species)
		herhoStandAlone(herhoTally.out)
}

