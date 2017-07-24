#ifndef ENCODING_HPP_INCLUDED__
#define ENCODING_HPP_INCLUDED__

#ifdef FFMPEG

extern "C" {
#include <libavutil/avassert.h>
#include <libavutil/channel_layout.h>
#include <libavutil/opt.h>
#include <libavutil/mathematics.h>
#include <libavutil/timestamp.h>
#include <libavutil/imgutils.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
#include <libswresample/swresample.h>
}

constexpr int BITRATE = 50000000;

class Mpeg {
  private:
    bool is_started =false;
    int framerate = 25;
    AVOutputFormat *fmt;
    AVFormatContext *oc;

    AVStream *st;
    AVCodecContext *enc;

    AVFrame *frame;
    AVFrame *rgb_frame;
    uint8_t *rgb_buffer;

    struct SwsContext *sws_ctx;

    AVFrame *alloc_picture(enum AVPixelFormat pix_fmt)
    {
        AVFrame *picture;

        picture = av_frame_alloc();
        if (!picture)
            return NULL;

        picture->format = pix_fmt;
        picture->width  = width;
        picture->height = height;

        av_frame_get_buffer(picture, 32);
        return picture;
    }

  public:
    int width;
    int height;
    bool IsStarted() { return is_started; }
    int AddFrame() {
        int ret;
        int got_packet = 0;
        AVPacket pkt = { 0 };

        glReadPixels (0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rgb_buffer);
        av_image_fill_arrays(rgb_frame->data, rgb_frame->linesize, rgb_buffer, AV_PIX_FMT_RGB24, width, height, 1);


        if (av_frame_make_writable(frame) < 0)
            return 1;

        // The picture is upside down - flip it:
        auto data = rgb_frame->data[0] + 3*width*height;
        uint8_t *flipped_data[4] = { data, data, data, data };
        int flipped_stride[4] = { -rgb_frame->linesize[0], -rgb_frame->linesize[1], -rgb_frame->linesize[2], -rgb_frame->linesize[3] };
        sws_scale(sws_ctx, flipped_data, flipped_stride, 0, enc->height, frame->data, frame->linesize);


        av_init_packet(&pkt);

        got_packet = 0;
        ret = avcodec_send_frame(enc, frame);
        if (ret < 0)
        {
            cerr << "Error encoding video frame: " << endl;
            return(1);
        }

        ret = avcodec_receive_packet(enc, &pkt);
        if (!ret)
            got_packet = 1;
        if (ret == AVERROR(EAGAIN))
            return 0;

        if (ret < 0) {
            cerr << "Error encoding video frame: " << endl;
            return 1;
        }

        if (got_packet) {
            /* rescale output packet timestamp values from codec to stream timebase */
            av_packet_rescale_ts(&pkt, enc->time_base, st->time_base);
            pkt.stream_index = st->index;

            /* Write the compressed frame to the media file. */
            ret = av_interleaved_write_frame(oc, &pkt);
        } else {
            ret = 0;
        }

        if (ret < 0) {
            cerr << "Error while writing video frame: " << endl;
            return(1);
        }

        return 0;
    }

    int Start(string filename) {
        AVCodec *video_codec;
        if(is_started) {
            cerr << "Stream already started" << endl;
            return 1;
        }
        is_started = true;

        GLint dims[4] = {0};
        glGetIntegerv(GL_VIEWPORT, dims);
        width = dims[2];
        height= dims[3];

        width = int((width+1)/4)*4+4;
        height = 2 * (height/2);

        int ret;
        int i;

        av_register_all();

        avformat_alloc_output_context2(&oc, NULL, NULL, filename.c_str());
//         oc->preload= (int)(0.5*AV_TIME_BASE);
        oc->max_delay= (int)(0.7*AV_TIME_BASE);

        fmt = oc->oformat;

        if (fmt->video_codec != AV_CODEC_ID_NONE) {
            /* find the encoder */
            video_codec = avcodec_find_encoder(fmt->video_codec);
            if (!(video_codec)) {
                cerr << "Could not find encoder for '" << avcodec_get_name(fmt->video_codec) << "'" << endl;
                return 1;
            }

            st = avformat_new_stream(oc, NULL);
            if (!st) {
                cerr << "Could not allocate stream\n";
                return 1;
            }
            st->id = oc->nb_streams-1;
            enc = avcodec_alloc_context3(video_codec);
            if (!enc) {
                cerr << "Could not alloc an encoding context\n";
                return 1;
            }

            enc->codec_id = fmt->video_codec;

            enc->bit_rate = BITRATE;
            enc->width    = width;
            enc->height   = height;
            AVRational tb;
            tb.num=1;
            tb.den=framerate;
            st->time_base = tb;
            enc->time_base       = st->time_base;

            enc->gop_size      = 200;
            enc->pix_fmt       = AV_PIX_FMT_YUV420P;
            if (enc->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
                enc->max_b_frames = 3;
            }
            if (enc->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
                enc->mb_decision = 2;
            }

            if (oc->oformat->flags & AVFMT_GLOBALHEADER)
                enc->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

//             enc->flags |= CODEC_FLAG_QSCALE;
//             enc->global_quality = 1180;
        }
        else {
            cerr << "could not init codecs!" << endl;
            return 1;
        }

        AVDictionary *opt = NULL;

        /* open the codec */
        ret = avcodec_open2(enc, video_codec, &opt);
        av_dict_free(&opt);
        if (ret < 0) {
            cerr << "Could not open video codec" << endl;
            return 1;
        }

        /* allocate and init a re-usable frame */
        frame = alloc_picture(enc->pix_fmt);
        if (!frame) {
            cerr << "Could not allocate video frame\n";
            return 1;
        }

        /* copy the stream parameters to the muxer */
        ret = avcodec_parameters_from_context(st->codecpar, enc);
        if (ret < 0) {
            cerr << "Could not copy the stream parameters\n";
            return 1;
        }

        av_dump_format(oc, 0, filename.c_str(), 1);

        if (!(fmt->flags & AVFMT_NOFILE)) {
            ret = avio_open(&oc->pb, filename.c_str(), AVIO_FLAG_WRITE);
            if (ret < 0) {
                cerr << "Could not open " << filename << " : " << endl;
                return 1;
            }
        }

        ret = avformat_write_header(oc, &opt);
        if (ret < 0) {
            cerr << "Error occurred when opening output file: " << endl;;
            return 1;
        }

        rgb_frame = alloc_picture(AV_PIX_FMT_RGB24);
        rgb_buffer = new uint8_t[width*height*4];
        sws_ctx = sws_getContext( width, height, AV_PIX_FMT_RGB24,
                                 width, height, AV_PIX_FMT_YUV420P,
                                 SWS_BICUBIC, NULL, NULL, NULL );
    }

    void Stop() {
        av_write_trailer(oc);
        avcodec_free_context(&enc);
        av_frame_free(&frame);
        sws_freeContext(sws_ctx);
        if (!(fmt->flags & AVFMT_NOFILE))
            avio_closep(&oc->pb);
        avformat_free_context(oc);
        delete [] rgb_buffer;
        is_started = false;
    }
};
#endif // FFMPEG
#endif // ENCODING_HPP_INCLUDED__
