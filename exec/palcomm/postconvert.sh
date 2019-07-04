##for palgrav communication PNG are not accepted so I write all in tiff and then convert in png using:
#this scruipt is used for palgrave formating counstraint
for i in *.tiff ; do convert $i ${i%.*}.png ; done
for i in *.tiff ; do convert $i ${i%.*}.jpg ; done
mv hdr_check_count0_random.tiff 3_hdr_check_count0_random.tiff
mv hdr_check_count0_conformist.tiff 4_hdr_check_count0_conformist.tiff
mv hdr_check_count0_topfiveS.tiff 5_hdr_check_count0_topfiveS.tiff
mv hdr_check_count0_topfiveAC.tiff 6_hdr_check_count0_topfiveAC.tiff
mv hdr_check_count0_topfiveAC.png 6_hdr_check_count0_topfiveAC.png
mv hdr_check_count0_topfiveS.png 5_hdr_check_count0_topfiveS.png
mv hdr_check_count0_random.png 3_hdr_check_count0_random.png
mv hdr_check_count0_conformist.png 4_hdr_check_count0_conformist.png
mv hdr_check_count0_topfiveAC.jpg 6_hdr_check_count0_topfiveAC.jpg
mv hdr_check_count0_topfiveS.jpg 5_hdr_check_count0_topfiveS.jpg
mv hdr_check_count0_random.jpg 3_hdr_check_count0_random.jpg
mv hdr_check_count0_conformist.jpg 4_hdr_check_count0_conformist.jpg
