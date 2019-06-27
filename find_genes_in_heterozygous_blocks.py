import sys
import os

#Script usage python genes_in_heter_regions.py GFF_FILE.gff BED_HETEROZYG_REGIONS.bed BED_COVERAGE.bed VCF_FILE.vcf OUTPUT_VCF_FILE.vcf
# Script outputs the number of genes located in regions of bed file, and will generate vcf file of variants located in genes


##########################  This part parses gff file and brings to format
						#	gene2647 ['333210', '334454', 'NC_018296.1']
						#	gene2646 ['332291', '333073', 'NC_018296.1']
						#	gene2645 ['330821', '330988', 'NC_018296.1']
						#	gene2644 ['327627', '329486', 'NC_018296.1']
						#   Gene, its coordinates, and on which chromosome it is located



with open(sys.argv[1], "r+") as gff_file, open(sys.argv[2], "r+") as bed_heterozygous, open(sys.argv[3], "r+") as vcf_file, open(sys.argv[4],"w") as vcf_output:
	gff={}
	for line in gff_file:
		if not line.startswith("#"):
			line=line.rstrip().split("\t")
			if line[2]=="gene":
				ID=line[8].split(";")[0][3:]
				coordinates=[]
				coordinates.extend([line[3],line[4], line[0]])
				gff[ID]=coordinates
	#~ for k,v in gff.items():
		#~ print k,v
############################

############################   This part parses bed file with heterozygous regions and  bring to format

							#	('2367748', '2367752') NC_018292.1
							#	('738823', '739671') NC_018300.1
							#	('236860', '237674') NC_018302.1

							# coordinates of heterozygous regions, and the chromosome it is located in

	bed_het={}
	for line in bed_heterozygous:
			if not line.startswith("#"):
				line=line.rstrip().split("\t")
				if line[0] == "chromosome1":
					col1="NC_018292.1"
				elif line[0] == "chromosome2":
					col1="NC_018295.1"
				elif line[0]=="chromosome3":
					col1="NC_018296.1"
				elif line[0]=="chromosome4":
					col1="NC_018297.1"
				elif line[0]=="chromosome5":
					col1="NC_018298.1"
				elif line[0]=="chromosome6":
					col1="NC_018300.1"
				elif line[0]=="chromosome7":
					col1="NC_018301.1"
				elif line[0]=="chromosome8":
					col1="NC_018302.1"
				bed_het[line[1],line[2]]= col1
	#print len(bed_het)
	#~ for k,v in bed_het.items():
		#~ print k,v
#########################################


######################################### This part finds genes located in heterozygous regions
	phased_genes_file=open("phased_genes.txt","w")
	genes_in_heter_reg={}
	for coord_region,chrom in bed_het.iteritems():
		for  gene,coord_gene in gff.iteritems():
			if chrom == coord_gene[2] and int(coord_gene[0]) >= int(coord_region[0]) and int(coord_gene[1]) <= int(coord_region[1]):
				#~ print gene, coord_gene, coord_region, chrom, coord_gene[2]
				genes_in_heter_reg[gene]=coord_gene
				#print gene
				phased_genes_file.write("%s\nalt_%s\n"%(gene,gene))
	phased_genes_file.close()
	phased_genes=len(genes_in_heter_reg)
	print "Phased %s genes"%(phased_genes)

##################################


################################## This part parses vcf file

	vcf=[]
	for line in vcf_file:
		line=line.rstrip()
		#print line[0]
		if "chromosome1" in line:
			line=line.replace("chromosome1","NC_018292.1")
			vcf.append(line)
			#print "NC_018292.1"
		elif "chromosome2" in line:
			line=line.replace("chromosome2","NC_018295.1")
			vcf.append(line)
		elif "chromosome3" in line:
			line=line.replace("chromosome3","NC_018296.1")
			vcf.append(line)
		elif "chromosome4" in line:
			line=line.replace("chromosome4","NC_018297.1")
			vcf.append(line)
		elif "chromosome5" in line:
			line=line.replace("chromosome5","NC_018298.1")
			vcf.append(line)
		elif "chromosome6" in line:
			line=line.replace("chromosome6","NC_018300.1")
			vcf.append(line)
		elif "chromosome7" in line:
			line=line.replace("chromosome7","NC_018301.1")
			vcf.append(line)
		elif "chromosome8" in line:
			line=line.replace("chromosome8","NC_018302.1")
			vcf.append(line)
		else:
			vcf.append(line)
	#~ for el in vcf:
		#~ print el

############### Find variant in genes and generate vcf file
	snps_in_genes={}
	for line in vcf:
		if line.startswith("#"):
			vcf_output.write("%s\n"%(line))
		else:
			line=line.split("\t")
			for gene,coords in genes_in_heter_reg.iteritems():
				if coords[2]==line[0] and int(coords[0])<=int(line[1]) and int(coords[1])>=int(line[1]) and line[9].split(":")[0] == "0/1":
					line="\t".join(line)
					vcf_output.write("%s\n"%(line))
					if gene in snps_in_genes:
						snps_in_genes[gene].append(line[3]+"_"+line[4])
					else:
						snps_in_genes[gene]=[line[3]+"_"+line[4]]
	print "vcf file with snps in phased genes is written"


##################################################################################################
