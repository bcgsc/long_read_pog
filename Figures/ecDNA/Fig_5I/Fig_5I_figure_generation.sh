#-----------------------
#Creating Fig. 5I
#----------------------

#Make and activate circos conda environment:
conda create -n circos_env -c bioconda circos
conda activate circos_env

#Run circos_ecDNA.conf script (saves PNG and SVG images into Fig_5I_images folder):
circos -conf circos_ecDNA.conf -svg

#(Additional colouring, labels, and legend added post-image generation)

