# circos.conf

karyotype = circos.txt


chromosomes = Proteome; Metabolome; Phenotype
chromosomes_radius = Proteome:1r ; Metabolome: 1r ; Phenotype: 0.85r

<ideogram>

<spacing>
default = 0.005r
</spacing>

radius    = 0.85r


thickness = 10p
fill      = yes
stroke = vdgrey
max_col_num = 1000
max_row_num = 1000

# Label settings
show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default 
label_radius     = 1r + 10p
center = yes
label_size       = 30
label_parallel   = yes


</ideogram>
#Links between environmental and experimental metabolites
<links>
radius = 0.87r
crest  = 1
ribbon           = yes
twist             = yes
stroke_color     = dgrey
color            = vlred
# thickness = 10

bezier_radius        = 0r
bezier_radius_purity = 0.5

<link>
file          = circosLinks.txt
</link>
# 
</links>

# All of the plots including: text, tiles
<plots>
<plot>

file = tiles.txt

show = yes
type = tile


r1 = 0.98r
r0 = 0.78r


layers      = 1
margin      = 0.02u
orientation = in
layers_overflow=collapse

thickness   = 190
padding     = 1

stroke_thickness = 5
stroke_color     = white

</plot>


<plot>
file = text.txt

show = yes
color = black
type = text

r1 = 0.96r
r0 = 0.84r

label_size   = 20p
label_font   = condensed

padding  = 0p
rpadding = 6p
</plot>

<plot>
file = clusterText.txt

show = yes
color = black
type = text

r1 = 1.35r
r0 = 1.05r

label_size   = 10p
label_font   = condensed

padding  = 0p
rpadding = 0p
</plot>

<plot>
type = line
show = yes

file = lines.txt
r1 = 1.1
r0 = 1.04

# min = 0.99
# max = 1.01

color = vvdgrey
thickness = 10
</plot>

</plots>



## All of the line plots
# <plot>
# # water controls
# show  = yes
# type  = heatmap

# file  = waterLines.txt
# r1    = 1r + 50p
# r0    = 1r + 10p
# min = 0
# max = 2.999398
# # max   = 7.07
# # min   = 0


# color     = dblue
# thickness = 10

# </plot>

# <plot>
# # Healthy Corals
# show  = yes
# type  = heatmap
# r1    = 1r + 95p
# r0    = 1r + 55p
# file  = healthyLines.txt

# min = 0
# max = 2.999398
# # max   = 7.07
# # min   = 0

# # color     = dgreen
# thickness = 10

# </plot>

# <plot>
# # PB Corals
# show  = yes
# type  = heatmap

# file  = pbLines.txt
# # r1    = 1r + 200p
# r1    = 1r + 140p
# r0    = 1r + 100p

# min = 0
# max = 2.999398
# # max   = 7.07
# # min   = 0


# # color     = dyellow
# thickness = 10

# </plot>

# <plot>
# # Bleached Corals
# show  = yes
# type  = heatmap

# file  = bleachLines.txt
# r1    = 1r + 145p
# r0    = 1r + 185p
# min = 0
# max = 2.999398
# # max   = 7.07
# # min   = 0

# # color     = dred
# thickness = 10

# </plot>



# </plots>


# <backgrounds>
# <background>
# color = vvlgrey
# y1    = 1r + 200p
# y0    = 1r + 2p
# </background>
# </backgrounds>

################################################################
# The remaining content is standard and required. It is imported 
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
angle_offset*= -42
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>