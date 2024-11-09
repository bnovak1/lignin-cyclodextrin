# Combine representative configurations for clusters with the scatter plots showing the clusters.
# Reads the lignol abbreviation from the command line.
# 
# Usage:
#   ./process_cluster_config_images.sh <lignol_abbreviation>
#
# Arguments:
#   lignol_abbreviation: A string representing the abbreviation of the lignol compound.
#
# Example:
#   ./process_cluster_config_images.sh GG

# Read lignol abbreviation from command line
lignol=$1

# Base directory
BASE_DIR=../analysis/${lignol}/clustering

# Image number sequence
if [ $lignol == "GG" ]; then
    image_sequence=(2 3 4 5 6 7 8)
elif [ $lignol == "TGG" ]; then
    image_sequence=(1 2 3 4 5)
elif [ $lignol == "GG_BB" ]; then
    image_sequence=(1 2 3)
fi

###################################################################################################
# Trim all images
###################################################################################################

convert ${BASE_DIR}/clusters.png -trim ${BASE_DIR}/clusters.png

for i in ${image_sequence[@]}
do
    convert ${BASE_DIR}/cluster_configs_${i}.tga -trim ${BASE_DIR}/cluster_configs_${i}.png
done

###################################################################################################
# Scale configuration images
###################################################################################################

if [ $lignol == "GG" ]; then
    scale_sequence=(115% 100% 120% 100% 100% 100% 100%)
elif [ $lignol == "TGG" ]; then
    scale_sequence=(100% 100% 100% 100% 100%)
elif [ $lignol == "GG_BB" ]; then
    scale_sequence=(100% 100% 115%)
fi

for idx in ${!image_sequence[@]}
do
    i=${image_sequence[$idx]}
    convert ${BASE_DIR}/cluster_configs_${i}.png \
        -resize ${scale_sequence[$idx]} \
        ${BASE_DIR}/cluster_configs_${i}_scaled.png
done

###################################################################################################
# Add symbols next to images
###################################################################################################

# sequence of colors for symbols
if [ $lignol == "GG" ]; then
    color_sequence=(pink black orange red purple blue green)
elif [ $lignol == "TGG" ]; then
    color_sequence=(cyan blue green red purple)
elif [ $lignol == "GG_BB" ]; then
    color_sequence=(blue red green)
fi

for idx in ${!image_sequence[@]}
do
    i=${image_sequence[$idx]}
    base_image=${BASE_DIR}/cluster_configs_${i}_scaled.png
    symbol_image=${BASE_DIR}/symbols/circle_${color_sequence[$idx]}.png
    convert $base_image $symbol_image \
        -gravity NorthWest \
        -composite ${BASE_DIR}/cluster_configs_${i}_symbol.png
done

###################################################################################################
# Combine images
###################################################################################################

border_width=20

if [ $lignol == "GG" ]; then

    # 5, 6 combined horizontally
    convert ${BASE_DIR}/cluster_configs_5_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_6_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/cluster_configs_5_6.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_5_6.png \
            -trim ${BASE_DIR}/cluster_configs_5_6.png

    # 7, 8, 4 combined vertically
    convert ${BASE_DIR}/cluster_configs_7_symbol.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/cluster_configs_8_symbol.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/cluster_configs_4_symbol.png \
            -bordercolor white -border 0x${border_width} \
            -append ${BASE_DIR}/cluster_configs_7_8_4.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_7_8_4.png \
            -trim ${BASE_DIR}/cluster_configs_7_8_4.png
    
    # 2, 3 combined horizontally
    convert ${BASE_DIR}/cluster_configs_2_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_3_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/cluster_configs_2_3.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_2_3.png \
            -trim ${BASE_DIR}/cluster_configs_2_3.png

    # 5-6, scatter plot, 2-3 combined vertically
    convert ${BASE_DIR}/cluster_configs_5_6.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/clusters.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/cluster_configs_2_3.png \
            -bordercolor white -border 0x${border_width} \
            -append ${BASE_DIR}/cluster_configs_5_6_scatter_2_3.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_5_6_scatter_2_3.png \
            -trim ${BASE_DIR}/cluster_configs_5_6_scatter_2_3.png

    # 7-8-4, 5-6-scatter-2-3 combined horizontally
    convert ${BASE_DIR}/cluster_configs_7_8_4.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_5_6_scatter_2_3.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/cluster_configs_7_8_4_5_6_scatter_2_3.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_7_8_4_5_6_scatter_2_3.png \
            -trim ${BASE_DIR}/cluster_configs_7_8_4_5_6_scatter_2_3.png

    # Rename the final image to clusters_configs.png
    mv ${BASE_DIR}/cluster_configs_7_8_4_5_6_scatter_2_3.png ${BASE_DIR}/clusters_configs.png

    # Remove temporary files
    rm ${BASE_DIR}/cluster_configs_5_6.png \
        ${BASE_DIR}/cluster_configs_7_8_4.png \
        ${BASE_DIR}/cluster_configs_2_3.png \
        ${BASE_DIR}/cluster_configs_5_6_scatter_2_3.png

elif [ $lignol == "TGG" ]; then

    # 5, 1 combined horizontally, add border for alignment
    convert ${BASE_DIR}/cluster_configs_5_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_1_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/cluster_configs_5_1.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_5_1.png \
            -trim ${BASE_DIR}/cluster_configs_5_1.png

    # 4, 3, 2 combined vertically
    convert ${BASE_DIR}/cluster_configs_4_symbol.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/cluster_configs_3_symbol.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/cluster_configs_2_symbol.png \
            -bordercolor white -border 0x${border_width} \
            -append ${BASE_DIR}/cluster_configs_4_3_2.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_4_3_2.png \
            -trim ${BASE_DIR}/cluster_configs_4_3_2.png

    # 5-1, scatter plot vertically
    convert ${BASE_DIR}/cluster_configs_5_1.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/clusters.png \
            -bordercolor white -border 0x${border_width} \
            -append ${BASE_DIR}/cluster_configs_5_1_scatter.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_5_1_scatter.png \
            -trim ${BASE_DIR}/cluster_configs_5_1_scatter.png

    # 4-3-2, 5-1-scatter combined horizontally
    convert ${BASE_DIR}/cluster_configs_4_3_2.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_5_1_scatter.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/clusters_configs.png

    # Trim the combined image
    convert ${BASE_DIR}/clusters_configs.png \
            -trim ${BASE_DIR}/clusters_configs.png

    # Remove temporary files
    rm ${BASE_DIR}/cluster_configs_5_1.png \
       ${BASE_DIR}/cluster_configs_4_3_2.png \
       ${BASE_DIR}/cluster_configs_5_1_scatter.png

elif [ $lignol == "GG_BB" ]; then

    # 1, 3, 2 combined horizontally
    convert ${BASE_DIR}/cluster_configs_1_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_3_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            ${BASE_DIR}/cluster_configs_2_symbol.png \
            -bordercolor white -border ${border_width}x0 \
            +append ${BASE_DIR}/cluster_configs_1_3_2.png

    # Trim the combined image
    convert ${BASE_DIR}/cluster_configs_1_3_2.png \
            -trim ${BASE_DIR}/cluster_configs_1_3_2.png

    # 1-3-2, scatter combined vertically
    convert ${BASE_DIR}/cluster_configs_1_3_2.png \
            -bordercolor white -border 0x${border_width} \
            ${BASE_DIR}/clusters.png \
            -bordercolor white -border 0x${border_width} \
            -append ${BASE_DIR}/clusters_configs.png

    # Trim the combined image
    convert ${BASE_DIR}/clusters_configs.png \
            -trim ${BASE_DIR}/clusters_configs.png

    # Remove temporary files
    rm ${BASE_DIR}/cluster_configs_1_3_2.png

fi