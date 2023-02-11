#!/bin/sh

# Create the star trail animation from the individual png images.

# Steps taken in the ffmpeg lines below (first block including title and explanation frames):
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
#
# Second block does not include the title frames and makes a video suitable for smooth looping.

USAGE="Usage: makevideo [-k] [-h]"
USAGELONG="Usage: makevideo [-k] [-a] [-t] [-h]\n Produce star trails video, by default suitable for smooth looping\n -k Produce 4K video\n -a Use HAP video codec\n -t Include title and explanation frames. Resulting video is not suitable for smooth looping\n -h help\n"
RESOLUTION="1920x1080"
IMFOLDER="images-2k"
CODEC="libx264"
PIX_FMT="yuv420p"
EXTENSION="mp4"
COMPR=""
TITLES=false

while getopts "kacth" options;
do
    case $options in
        k) RESOLUTION="3840x2160"
            RESFACTOR=2
            IMFOLDER="images-4k"
            ;;
        a) CODEC="hap"
            PIX_FMT="rgba"
            EXTENSION="mov"
            COMPR="-compressor snappy"
            ;;
        t) TITLES=true
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

if $TITLES
then
    ffmpeg \
        -f lavfi -i "color=c=black:s=${RESOLUTION}:r=25" \
        -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/edr3-title-frame.png" \
        -loop 1 -framerate 25 -t 7 -i "${IMFOLDER}/star-trail-text-frame.png" \
        -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/startframe.png" \
        -framerate 25 -i "${IMFOLDER}/frame%04d.png" \
        -loop 1 -framerate 25 -t 3 -i "${IMFOLDER}/endframe.png" \
        -filter_complex \
        "[1:v]fade=type=out:duration=1:start_time=2,format=${PIX_FMT}[v0]; \
        [2:v]fade=type=in:duration=1,fade=type=out:duration=1:start_time=6,format=${PIX_FMT}[v1]; \
        [0:v][3:v]overlay=shortest=1,format=${PIX_FMT}[v2]; \
        [0:v][4:v]overlay=shortest=1,tpad=stop_mode=clone:stop_duration=3,format=${PIX_FMT}[v3]; \
        [0:v][5:v]overlay=shortest=1,fade=type=in:duration=1,format=${PIX_FMT}[v4]; \
        [v0][v1][v2][v3][v4]concat=n=5" -vcodec $CODEC -s $RESOLUTION $COMPR trails.$EXTENSION
else
    ffmpeg \
        -f lavfi -i "color=c=black:s=${RESOLUTION}:r=25" \
        -loop 1 -framerate 25 -t 2 -i "${IMFOLDER}/panoramaframe.png" \
        -loop 1 -framerate 25 -t 4 -i "${IMFOLDER}/startframe.png" \
        -framerate 25 -i "${IMFOLDER}/frame%04d.png" \
        -loop 1 -framerate 25 -t 6 -i "${IMFOLDER}/endframe.png" \
        -loop 1 -framerate 25 -t 2 -i "${IMFOLDER}/panoramaframe.png" \
        -filter_complex \
        "[0:v][1:v]overlay=shortest=1,format=${PIX_FMT}[v0]; \
        [0:v][2:v]overlay=shortest=1,fade=type=in:duration=1,format=${PIX_FMT}[v1]; \
        [0:v][3:v]overlay=shortest=1,tpad=stop_mode=clone:stop_duration=3,format=${PIX_FMT}[v2]; \
        [0:v][4:v]overlay=shortest=1,fade=type=in:duration=1,fade=type=out:start_time=5:duration=1,format=${PIX_FMT}[v3]; \
        [0:v][5:v]overlay=shortest=1,fade=type=in:duration=1,format=${PIX_FMT}[v4]; \
        [v0][v1][v2][v3][v4]concat=n=5" -vcodec $CODEC -s $RESOLUTION $COMPR trails.$EXTENSION
fi
