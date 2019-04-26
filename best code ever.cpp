BYTE * ikinciproje(BYTE *image, int w, int h) {
	int  h1 = h + 2;
	int  w1 = w + 2;
	int hold=0;
	BYTE * buffer2 = new BYTE[w1*h1];
	BYTE * buffer3 = new BYTE[w1*h1];
	for (int i = 0; i < w1*h1; i++)
		{			
		buffer2[i] = 0;
		buffer3[i] = 0;
		}
	
	for (int j = 1; j < h1 - 1; j++) {
		for (int i = 1; i < w1 - 1; i++) {
			buffer2[j*w1 + i] = image[hold];
			buffer3[j*w1 + i] = image[hold];
			hold++;
		}
	}
	
	BYTE * ximage = new BYTE[w1*h1];
	
	BYTE * yimage = new BYTE[w1*h1];
	
	BYTE * totalimage = new BYTE[w1*h1];
	
	
	double total = 0;
	double total2=0;
	double mask1[9] = { 1,2,1,0,0,0,-1,-2,-1 };
	double mask2[9] = { 1,0,-1,2,0,-2,1,0,-1 };
	for (int raw = 1; raw < h1 -1; raw++) {
		for (int col = 1; col < w1 - 1; col++) {

			for (int i = -1; i < 2; i++) {
				for (int j = -1; j <2; j++) {
					total += buffer2[(raw + i)*w1 + (col + j)] * mask1[(i +1) * 3 + j + 1];
					total2 += buffer3[(raw + i)*w1 + (col + j)] * mask2[(i + 1) * 3 + j + 1];
				}
			}
			ximage[raw*w1 + col] =(int)abs(total);
			yimage[raw*w1 + col] = (int)abs(total2);
			total = 0;
			total2 = 0;
		}
	}
	
	BYTE * holdangle = new BYTE[w1*h1];
	int aci = 0;

	
	for (int i = 0; i < w1*h1; i++) {
	
			totalimage[i] =ximage[i]+yimage[i];
			if (ximage[i] == 0 && yimage[i] == 0)
				holdangle[i] = 0;
			else if (ximage[i] == 0 && yimage[i] != 0)
				holdangle[i] = 90;
			else {
				aci = atan(yimage[i]/ximage[i]  );
				if (((aci < 22.5) && (aci > -22.5)) || (aci > 157.5) || (aci < -157.5))
					holdangle[i] = 0;
				else if (((aci > 22.5) && (aci < 67.5)) || ((aci < -112.5) && (aci > -157.5)))
					holdangle[i] = 45;
				else if (((aci > 67.5) && (aci < 112.5)) || ((aci < -67.5) && (aci > -112.5)))
					holdangle[i] = 90;
				else if (((aci > 112.5) && (aci < 157.5)) || ((aci < -22.5) && (aci > -67.5)))
					holdangle[i] = 135;

			
			}

			if (totalimage[i] < 150)
				totalimage[i] = 0;
			else totalimage[i] == 255;
		}
	
	return totalimage;
}