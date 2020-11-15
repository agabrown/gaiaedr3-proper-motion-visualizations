#!/bin/sh

# Create the star trail animation from the individual png images.

# Steps taken in the ffmpeg lines below
#
# 1 ffmpeg command
# 2 generate black backgrond image (stream [0:v])
# 3 title frame is kept static for 3 seconds (stream [1:v])
# 4 frame with text explanation is kept static for 7 seconds (stream [2:v])
# 5 start frame with star positions is static for 3 seconds (stream [3:v])
# 6 input frames for full animation (stream [4:v])
# 7 start complex filter operating on the above streams
# 8 fade out title frame during 1 second (starting after 2 seconds) and store in [v0]
# 9 fade in text frame and fade it out after 6 seconds, store in [v1]
# 10 put black background behind the start frame, store in [v2]
# 11 put black background behind the animation frames and keep last frame frozen for 3 seconds, store in [v3]
# 12 put black background behind the endframe, store in [v4]
# 13 concatenate the streams v0-v4

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
