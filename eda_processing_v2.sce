clf() // reiniciar graficos


// EXTRAIR DADOS E DEFINIR SINAL
sinal=csvRead('eda_raw_pd.csv')
L=length(sinal)
fs=1000 // samp. freq 1kHz
k1=217 // size of filter 1
k2=6000 // size of filter 2
samples = 1:L
time = (1/1000)*samples

subplot(3,2,1)
title("Sinal bruto",'fontsize',8)
plot(time,sinal,'b') // sinal original, com ruido


// OBTER CAUSA DO RUIDO => conclui-se 50 Hz com harmonicas nos n50Hz
S=fft(sinal-mean(sinal))
dnw=2*%pi/length(S)
faxis=(fs/2/%pi)*[0:dnw:2*%pi-dnw]

subplot(3,2,2)
title("FFT do sinal",'fontsize',6)
plot(S,'r')
subplot(3,2,3)
title("Harmónicas do sinal (50Hz)",'fontsize',7)
plot(faxis,abs(S))


// APLICAR PRIMEIRO FILTRO LP PARA REMOVER RUIDO DO SINAL ORIGINAL
before = repmat(sinal(1),k1,1)
after = repmat(sinal($),k1,1) // $ - last index
sinalr = cat(1,before,sinal,after)
h1=(1/k1)*ones(1,k1)
w1=filter(h1,1,sinalr)
w1(1:(1.5*k1))=[]
w1(L+1:L+1+k1)=[]

subplot(3,2,4)
title("Filtro passa-baixo", "fontsize", 8)
plot(time,w1,"r",'LineWidth',4)


// APLICAR SEGUNDO FILTRO LP PARA EXTRAIR SCL
before = repmat(w1(1),k2,1)
after = repmat(w1($),k2,1)
w1r = cat(1,before,w1,after)
h2=(1/k2)*ones(1,k2)
w2=filter(h2,1,w1r)
w2(1:(1.5*k2))=[]
w2(L+1:L+1+k2)=[]

subplot(3,2,5)
title("Skin Conductance Level", "fontsize", 8)
plot(time,w1,"r",'LineWidth',3)
plot(time,w2,"b",'LineWidth',3)


// EXTRAIR SCR
dw=w1-w2
scr=dw(2000:length(dw))
scr=dw

subplot(3,2,6)
title("Skin Conductance Response", "fontsize", 8)
plot(time,scr,'LineWidth',2)


// SAVE DATA
write('eda_filtered_scilab.csv',format(20),w1)
write('scl_scilab.csv',format(20),w2)
write('scr_scilab.csv',format(20),scr)

