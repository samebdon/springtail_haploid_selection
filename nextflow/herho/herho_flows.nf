include {herhoTally; herhoStandAlone} from './herho_tasks.nf'

workflow herho_flow {
	take:
		vcf
		bed
		species
	main:
		herhoTally(vcf, bed, species)
		herhoStandAlone(herhoTally.out)
}

