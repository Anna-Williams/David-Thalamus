library(pals) # for palettes with large n #kelly()22, #polychrome()#36, cols

# remove the black and white from the pallete, and the light blue similars to cols25
# assessed with pal.bands and pal.cube
kelly_col <- unname(kelly()[-c(1,2,6)])
# remove the orange colour that is similar to the kelly, and sort with prettiest colours first and ugliest or still similar to kelly last
cols24 <- unname(cols25()[c(19,22:24,8:14, 1:4, 15:17,6,7,25,18,20,21)])
# merge all pallettes for long list colours, 
# my favourite are the first 22 (kelly), cols24 is not bad and I tried to keep it distinct
# and polychrome is just added at the end in case we are missing levels (no checked it is safe with the other 42)
cols <- c(kelly_col, cols24, unname(polychrome()))
# These colours are taken from D. Munro images, and were adjusted with pal.safe
col_wt_het_ko <- c("#6c6c6c", "#7c83b9", "#e28c93")

# If we need 6 groups, playing with muted, still safe ( limit between HET and WT in 
#BW, but pretty decent)
col_2x_wt_het_ko <- c("#444444",scales::muted(col_wt_het_ko)[2], scales::muted(col_wt_het_ko[3], l=40),  col_wt_het_ko)
# These colours were chosen by me with pal.safe at first
col_wt_ko <- c("#666666", "#E25822")

# #at the end we decided to keep consistent with young ones
col_wt_ko <- c("#6c6c6c", "#e28c93")