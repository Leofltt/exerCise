/*

	biquad.c

	an implementation of a biquad filter using portaudio
	+ 
	some fun feedback action just because	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "portaudio.h" 
#define SR 44100
#define PI 3.1415926535897
#define BUF 256

//STRUCT FOR FILTER DATA
typedef struct {
	float a0, a1, a2, b0, b1, b2, x_uno, x_due, y_uno, y_due;
	char* type;
}biquad;

//STRUCT FOR FEEDBACK
typedef struct { 
	float delay[(int)SR/2];
	int rp;
	float fb;
	 
}echo;


biquad* new_filt(int filt_type, float freq, float Q, float gain, int samp_rate);

float compute_filt(float input, biquad* bq);
void coef_filt(biquad *bq, int filt_type, float omega, float sn, float c, float A, float alpha, float beta);
void free_filt(biquad* bq);
int pa_callBack(const void *input, void *output, unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* time_info,
				PaStreamCallbackFlags flags,
				void *userData);

biquad *new_filt(int filt_type, float freq, float Q, float gain, int samp_rate) {
	
	biquad *temp= (biquad *)malloc(sizeof(biquad));

	float A = pow(10, gain / 40); 
    float omega = 2 * PI * freq / samp_rate;
    float sn = sin(omega);
    float c = cos(omega);
    float alpha = sn / (2*Q);
    float beta = sqrt(A + A);

    coef_filt(temp, filt_type, A, omega, sn, c, alpha, beta);
    
    //scale by a0
    temp->a1 = temp->a1 / temp->a0;
	temp->a2 = temp->a2 / temp->a0;
	temp->b0 = temp->b0 / temp->a0;
	temp->b1 = temp->b1 / temp->a0;
	temp->b2 = temp->b2 / temp->a0;

	//init history
	temp->x_uno = 0.0;
	temp->x_due = 0.0;
	temp->y_uno = 0.0;
	temp->y_due = 0.0;

	return temp;
}

void coef_filt(biquad *bq, int filt_type, float omega, float sn, float c, float A, float alpha, float beta) {
	
	switch (filt_type) {
		
		//LOWPASS COEFFICIENTS
		case 1:
			bq->type = "lpf";
    		bq->a0 = 1.0 + alpha;
    		bq->a1 = -2.0 * c;
    		bq->a2 = 1.0 - alpha;
    		bq->b0 = (1.0 - c) /2.0;
    		bq->b1 = 1.0 - c;
    		bq->b2 = (1.0 - c) /2.0;
			break;

		//HIGHPASS COEFFICIENTS	
		case 2:
			bq->type = "hpf";
    		bq->a0 = 1 + alpha;
    		bq->a1 = -2 * c;
    		bq->a2 = 1 - alpha;
    		bq->b0 = (1 + c) /2.0;
    		bq->b1 = -(1 + c);
    		bq->b2 = (1 + c) /2.0;
			break;

		//BANDPASS COEFFICIENTS
		case 3:
			bq->type = "bpf";
			bq->b0 = alpha;
    		bq->b1 = 0;
    		bq->b2 = -alpha;
    		bq->a0 = 1 + alpha;
    		bq->a1 = -2 * c;
    		bq->a2 = 1 - alpha;
			break;

		//LOWSHELF COEFFICIENTS
		case 4:
			bq->type = "lsf";
			bq->b0 = A * ((A + 1) - (A - 1) * c + beta * sn);
    		bq->b1 = 2 * A * ((A - 1) - (A + 1) * c);
    		bq->b2 = A * ((A + 1) - (A - 1) * c - beta * sn);
    		bq->a0 = (A + 1) + (A - 1) * c + beta * sn;
    		bq->a1 = -2 * ((A - 1) + (A + 1) * c);
    		bq->a2 = (A + 1) + (A - 1) * c - beta * sn;
			break;
		
		//HIGHSHELF COEFFICIENTS
		case 5:
			bq->type = "hsf";
			bq->a0 = (A + 1) - (A - 1) * c + beta * sn;
    		bq->a1 = 2 * ((A - 1) - (A + 1) * c);
    		bq->a2 = (A + 1) - (A - 1) * c - beta * sn;
			bq->b0 = A * ((A + 1) + (A - 1) * c + beta * sn);
    		bq->b1 = -2 * A * ((A - 1) + (A + 1) * c);
    		bq->b2 = A * ((A + 1) + (A - 1) * c - beta * sn);
			break;	

		//PEAK COEFFICIENTS	
		case 6:
			bq->type = "peak";
			bq->a0 = 1 + (alpha /A);
    		bq->a1 = -2 * c;
    		bq->a2 = 1 - (alpha /A);
			bq->b0 = 1 + (alpha * A);
    		bq->b1 = -2 * c;
    		bq->b2 = 1 - (alpha * A);
			break;

		//NOTCH COEFFICIENTS	
		case 7:
			bq->type = "notch";
			bq->a0 = 1 + alpha;
    		bq->a1 = -2 * c;
    		bq->a2 = 1 - alpha;
			bq->b0 = 1;
    		bq->b1 = -2 * c;
    		bq->b2 = 1;
			break;				
	}
}

float compute_filt(float input, biquad* bq) {
	
	//compute  output
	float output = 	(bq->b0 * input) + (bq->b1 * bq->x_uno) + (bq->b2* bq->x_due) - (bq->a1* bq->y_uno) - (bq->a2* bq->y_due);
	
	//update sample history
	bq->x_due = bq->x_uno;
	bq->x_uno = input;
	bq->y_due = bq->y_uno;
	bq->y_uno = output;
	return output;
}

void free_filt(biquad* bq) {
	free(bq);
}

biquad *bq;
int main(int argc, char** argv){
	
	PaError err;
	err = Pa_Initialize();

	//handle
	PaStream *stream;
	echo *data = (echo *) calloc(sizeof(echo),1);
	
	int ft;
	float fr, qu, gaindb;

	if (argc != 6) {
		fprintf(stderr, "\nBIQUAD FILTERING OF MIC INPUT \nUSAGE: ./name filter type, frequency, Q, gain (db), feedback \n  \nfilter types: \n 1 = lowpass \n 2 = highpass \n 3 = bandpass \n 4 = lowshelf \n 5 = highshelf \n 6 = peak \n 7 = notch\n" );
		return 1;
	} else {
		ft = atoi(argv[1]);
		fr = atof(argv[2]);
		qu = atof(argv[3]);
		gaindb = atof(argv[4]);
		data->fb = atof(argv[5]);
	}

	bq = new_filt(ft,fr,qu,gaindb,SR);
	err = Pa_OpenDefaultStream(	&stream,
								1,
								1,
								paFloat32,
								SR,
								BUF,
								pa_callBack,
								data);

	err = Pa_StartStream(stream);

	while(1);

	err = Pa_StopStream(stream);
	err = Pa_CloseStream(stream);
	err = Pa_Terminate();

	free_filt(bq);
	free(data);
	return 0;
}

int pa_callBack(const void *input, void *output, unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* time_info,
				PaStreamCallbackFlags flags,
				void *userData) {
	echo *data = (echo *) userData;
	float out;
	float *inp = (float *)input, *outp = (float *)output, *delay = data->delay;
	int i;
	int rp = data->rp;
	for(i = 0; i < framesPerBuffer; i++)
	{
		out = delay[rp];
		delay[rp++] = compute_filt(*inp++, bq) + out * data->fb;
		if(rp >= SR/2) rp = 0;
		*outp++ = out;
	}
	data->rp = rp;
	return 0;
}
