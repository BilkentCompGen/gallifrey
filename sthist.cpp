#include "sthist.h"
#include <iostream>

#include "histogram.h"

#define TB_IMPL
#define TB_OPT_V1_COMPAT
#include "ext/termbox2.h"
int run_hist_on_stdin(int dims, int freq = 10) {
    int w = 1;
    int h = 1;

    while (w * h < dims) {
        if (w > h)
            ++h;
        else if (w <= h)
            ++w;
    };
    Screen screen {h,w};
    w = 0;
    h = 0;
    int histograms [dims];
    int idx = 0;
    while (w * h < dims) {
        histograms[idx++] = screen.addHistogram<double>(h,w, "h");
        if (w > h){
            ++h;
        }   
        else if (w <= h){
            ++w;
        }
            
    };
    
    double val;
    int i = 0;
    struct tb_event ev;

    while(std::cin >> val){
        screen.pushHistogram(histograms[i%dims], val);
        ++i;
        if(i%freq==0){
            screen.update();
        }
        while (tb_peek_event(&ev, 0) != TB_ERR_NO_EVENT) {
            if (ev.type == TB_EVENT_KEY && (ev.key == TB_KEY_CTRL_C || ev.key == TB_KEY_ESC || ev.key == TB_KEY_ENTER))
                goto enough;
        }
    }
    enough:
    while (tb_peek_event(&ev, 100) == TB_ERR_NO_EVENT) {
        if (ev.type == TB_EVENT_KEY)
            goto enough;
    }
    return 0;
}