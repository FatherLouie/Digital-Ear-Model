# Digital-Ear-Model
This project implements the auditory model of the human ear as series of several discrete time filters. The premise of this exercise comes from the course project
of the DSP course at UNSW Sydney, the extensive details of which can be found here: [Cochlear Signal Processing](https://ieeexplore.ieee.org/document/9054297).
<br><br>
The MATLAB code presented here implements said digital model. For the outer and middle ears, the magnitude of the frequency response is well known experimentally. The filters have been created in order to best fir this experimental data.
<br><br>
The cochlea (inner ear) has been implemented as a series of 128 cascading filters. It has been verified by checking the excitation of different regions of the cochlea on passing various input signals. All plots for the mentioned tests are included in the MATLAB code that has been written.
