#!/bin/sh

# Create the star trail animation from the individual png images.

# Steps taken in the ffmpeg lines below
#
# 1 ffmpeg command
# 2 generate black backgrond image (stream [0:v])
# 3 first frame is static for a few seconds (stream [1:v])
# 4 input frames for full animation (stream [2:v])
# 5 start complex filter operating on the three above streams
# 6 put background behind static frame (stream [v0])
# 7 put background behind animation (stream [v1])
# 8 concatenate the two previous streams

USAGE="Usage: makevideo [-k] [-h]"
USAGELONG="Usage: makevideo [-k] [-h]\n -k produce 4K video\n -h help\n"
RESOLUTION="1920x1080"
IMFOLDER="images-2k"

while getopts "kh" options;
do
    case $options in
        k) RESOLUTION="3840x2160"
            IMFOLDER="images-4k"
            ;;
        h)
            echo -e $USAGELONG
            exit 0
            ;;
        \?)
            echo $USAGE
            exit 1
            ;;
    esac
done
shift $(($OPTIND-1))


ffmpeg \
    -f lavfi -i "color=c=black:s=${RESOLUTION}:r=25" \
    -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/edr3-title-frame.png" \
    -loop 1 -framerate 25 -t 7 -i "${IMFOLDER}/star-trail-text-frame.png" \
    -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/startframe.png" \
    -framerate 25 -i "${IMFOLDER}/frame%04d.png" \
    -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/endframe.png" \
    -filter_complex \
    "[1:v]fade=type=out:duration=1:start_time=2,format=yuv420p[v0]; \
    [2:v]fade=type=in:duration=1,fade=type=out:duration=1:start_time=6,format=yuv420p[v1]; \
    [0:v][3:v]overlay=shortest=1,format=yuv420p[v2]; \
    [0:v][4:v]overlay=shortest=1,tpad=stop_mode=clone:stop_duration=3,format=yuv420p[v3]; \
    [0:v][5:v]overlay=shortest=1,fade=type=in:duration=1,format=yuv420p[v4]; \
    [v0][v1][v2][v3][v4]concat=n=5" -vcodec libx264 -s $RESOLUTION trails.mp4
