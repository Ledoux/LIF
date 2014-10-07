###LOAD PACKAGES###
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
import itertools
import scipy


####C SWIG WRAPPING########
def getFilterDictByType(**Dict):
    FilteredDict={'DoubleDict':{},'IntDict':{},'StringDict':{}};
    for Key,Value in Dict.items():
        TypeDictName=getCTypeNameFromPythonType(type(Value))+'Dict';
        FilteredDict[TypeDictName][Key]=Value;
    return FilteredDict;
def getCTypeNameFromPythonType(PythonType):
    if PythonType in [float]:
        return 'Double';
    elif PythonType in [int]:
        return 'Int';
    elif PythonType in [str]:
        return 'String'; 
def getCArgsFromDict(Dict):
    CArgs=[]
    DictOrderedKeys=Dict.keys()
    for Key in sorted(Dict):
        CArgs.append(Dict[Key]);
    return CArgs;


#####GENERIC PYTOOLCLASS
class PyToolClass():

	#########################################################################
	########GET ACCESS TO ATTRIBUTES LIKE A DICTIONNARY/JSON OBJECT#########
	def __getitem__(self, key): return self.__dict__[key];
	def __setitem__(self, key, item):self.__dict__[key]=item;
	def __delitem__(self,key):del self.__dict__[key];
	
	def keys(self):
		return self.__dict__.keys();
	def values(self):
		return self.__dict__.values();
	def items(self):
		return self.__dict__.items();
		
	def printAll(self):
		print(self.__dict__);
        
	def printParam(self,Name):
		print("%s :"%Name);
		print(self[Name]);


####COMPLEX NUMBER METHODS

def getArg(z):
	return 2.*np.arctan(np.imag(z)/(np.sqrt(np.imag(z)**2+np.real(z)**2)+np.real(z)));

def combi(k,n):
	if(k<=n):
		return math.factorial(n)/(math.factorial(n-k)*math.factorial(k));
	else:
		return 0;

def getCombi(stuff):
	subsubsubset=list();
	for L in range(0, len(stuff)+1):
		subsubset=list();
		for subset in itertools.combinations(stuff, L):
			subsubset.append(subset);
		subsubsubset.append(subsubset);
	return subsubsubset;

def write2D(data,fname):
	output = open(fname, 'wb');
	dim=shape(data);
	j=0;
	for z in range(0,dim[0]):
		k=0;
		for y in range(0,dim[1]):
			output.write('%.3f ' %data[j,k]);
			k=k+1;	
		output.write('\n');	
		j=j+1;	
	output.close();

def read2D(fname):
	a=open(fname,'r');
	data=list();
	li=a.readline().rsplit(' ');
	while len(li) > 1 :
		data.append(array(li)[array(range(0,len(li)-1))]);
		li=a.readline().rsplit(' ');
	return array(data)

def getContour2(X,Y,ech):
	x0, y0 = np.mean(X), np.mean(Y);print(x0,y0);
	C = (X - x0) + 1j * (Y - y0);
	angles = np.angle(C);
	distances = np.absolute(C);
	sortidDist = [i[0] for i in sorted(enumerate(distances), key=lambda x:x[1],reverse=True)]
	anglesDist = angles[ sortidDist ][range(0,ech)];
	distancesDist = distances[ sortidDist ][range(0,ech)];
	XDist=X[ sortidDist ][range(0,ech)];
	YDist=Y[ sortidDist ][range(0,ech)];
	sortidAng = [i[0] for i in sorted(enumerate(anglesDist), key=lambda x:x[1],reverse=True)]
	anglesAng = anglesDist[ sortidAng];
	distancesAng = distancesDist[ sortidAng ];
	XAng=XDist[ sortidAng ];
	YAng=YDist[ sortidAng ];
	return XAng,YAng;

def getColor(nColors):
	jet = cm = plt.get_cmap('jet');
	cNorm  = colors.Normalize(vmin=0, vmax=nColors);
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet);
	return [scalarMap.to_rgba(i) for i in range(0,nColors)];

def getContour(X,Y,nEch):
	x0, y0 = np.mean(X), np.mean(Y);
	C = (X - x0) + 1j * (Y - y0);
	angles = arg(C);
	sortAng=[i[0] for i in sorted(enumerate(angles), key=lambda x:x[1])];
	anglesSort=angles[sortAng];
	CSort=C[sortAng];
	scanAng=arange(-pi,pi+2*pi/float(nEch),2*pi/float(nEch));
	CCont=list();
	#colorRange=getColor(len(scanAng));
	#colorRange=list();
	colorSet=list();
	CScan=[list() for i in range(0,len(scanAng))];
	who=[list() for i in range(0,len(scanAng))];
	whoMax=list();
	count=0;
	for i in range(0,len(scanAng)):
		ang=scanAng[i];
		#colorSet.append(colorRange[i]);
		if i%2==0:
			colorSet.append("black");
		else:
			colorSet.append("red");
		while count<len(anglesSort) and anglesSort[count]<ang:
			who[i].append(sortAng[count]);
			CScan[i].append(CSort[count]);
			count+=1;
		if len(CScan[i])>0:
			distances=np.absolute(array(CScan[i]));
			whereMax=where((distances)==max(distances))[0][0];
			whoMax.append(who[i][whereMax]);
			CCont.append(distances[whereMax]*np.exp(1j*np.angle(CScan[i][whereMax])));
	CCont.append(CCont[0]);
	CCont=array(CCont);

				
	return x0+real(CCont),y0+imag(CCont),[x0+real(CScan[i]) for i in range(0,len(scanAng))],[y0+imag(CScan[i]) for i in range(0,len(scanAng))],colorSet,whoMax;


def getContour3(X,Y,x0,y0,minTheta,maxTheta,nEch):
	C = (X - x0) + 1j * (Y - y0);
	angles = arg(C);
	sortAng=[i[0] for i in sorted(enumerate(angles), key=lambda x:x[1])];
	anglesSort=angles[sortAng];
	CSort=C[sortAng];
	deltaTheta=maxTheta-minTheta;
	scanAng=arange(minTheta,maxTheta+(deltaTheta/float(nEch)),deltaTheta/float(nEch));
	CCont=list();
	#colorRange=getColor(len(scanAng));
	#colorRange=list();
	colorSet=list();
	CScan=[list() for i in range(0,len(scanAng))];
	who=[list() for i in range(0,len(scanAng))];
	whoMax=list();
	count=0;
	for i in range(0,len(scanAng)):
		ang=scanAng[i];
		#colorSet.append(colorRange[i]);
		if i%2==0:
			colorSet.append("black");
		else:
			colorSet.append("red");
		while count<len(anglesSort) and anglesSort[count]<ang:
			who[i].append(sortAng[count]);
			CScan[i].append(CSort[count]);
			count+=1;
		if len(CScan[i])>0:
			distances=np.absolute(array(CScan[i]));
			whereMax=where((distances)==max(distances))[0][0];
			whoMax.append(who[i][whereMax]);
			CCont.append(distances[whereMax]*np.exp(1j*np.angle(CScan[i][whereMax])));
	CCont.append(CCont[0]);
	CCont=array(CCont);
	
	
	return x0+real(CCont),y0+imag(CCont),[x0+real(CScan[i]) for i in range(0,len(scanAng))],[y0+imag(CScan[i]) for i in range(0,len(scanAng))],colorSet,whoMax;


def getFirstNonZeroValues(number):
	decimal=0.;
	if number>0.:
		while number<1.:
			number=10.*number;
			decimal+=1;
		return number,decimal
	else:
		return 0.,-1



def PCA(data, dims_rescaled_data=2):
	""" performs principal components analysis
	(PCA) on the n-by-p data matrix A
	Rows of A correspond to observations, columns to variables. 

	Returns :
	coeff :
	is a p-by-p matrix, each column containing coefficients 
	for one principal component.
	score : 
	the principal component scores; that is, the representation 
	of A in the principal component space. Rows of SCORE 
	correspond to observations, columns to components.
	latent : 
	a vector containing the eigenvalues 
	of the covariance matrix of A.
	"""
	M = (data-mean(data.T,axis=1)).T # subtract the mean (along columns)
	[latent,coeff] = linalg.eig(cov(M)) # attention:not always sorted
	score = dot(coeff.T,M) # projection of the data in the new space
	return coeff,score,latent

def plot_pca(data,nameData,coeff,score,latent,pred=[],choose=[[0,1]],nEigMax=1000):
	figure();
	fontSize=25;
	plt.subplots(figsize=(50,50));
	nRows=35;
	nCols=30;
	colLen=4;
	rowLen=4;
	row=0;
	col=1;
	ax=list();
	ax.append(plt.subplot2grid((nRows,nCols), (row,col), rowspan=1,colspan=1));
	row+=2;
	ax=list();
	for j in xrange(len(choose)/2):
		ax.append(plt.subplot2grid((nRows,nCols), (row,col), rowspan=rowLen,colspan=colLen));
		col+=colLen+1;
	row+=rowLen+2;
	col=1;
	for j in xrange(len(choose)/2,len(choose)):
		ax.append(plt.subplot2grid((nRows,nCols), (row,col), rowspan=rowLen,colspan=colLen));
		col+=colLen+1;
	row+=rowLen+2;
	col=1;
	for j in xrange(2):
		ax.append(plt.subplot2grid((nRows,nCols), (row,col), rowspan=rowLen,colspan=colLen));
		col+=colLen+1;
	row+=rowLen+2;
	col=1;
	for j in xrange(1):
		ax.append(plt.subplot2grid((nRows,nCols), (row,col), rowspan=rowLen,colspan=colLen+5));
		col+=colLen+1;
		
	
	# every eigenvector describe the direction
	# of a principal component.
	# sorted them by eigenvalue in decreasing order
	idx = np.argsort(latent)[::-1];
	coeff = coeff[idx,:];
	#for i in xrange(len(idx)):
	#	coeff[i,:]/=sum(coeff[i,:]);
	latent = latent[idx];
	#nameData=nameData[idx];
	colorEig=[];
	colorEigTable=['red','orange','yellow','green','blue','violet','brown','black'];
	count=0;
	for i in xrange(len(idx)):
		colorEig.append(colorEigTable[count]);
		count+=1;
		if count==len(colorEigTable):
			count=0;

	m = mean(data,axis=1)
	lenD=100.;
	fig=-1;
	colorPoint="black";
	for i in xrange(len(choose)):
		fig+=1;
		ax[fig].plot(data[choose[i][0],:],data[choose[i][1],:],'o',color=colorPoint) # the data
		for k in xrange(min(len(latent),nEigMax)):
			ax[fig].plot([0,coeff[k,choose[i][0]]*lenD]+m[choose[i][0]], [0,coeff[k,choose[i][1]]*lenD]+m[choose[i][1]],'-',lw=3,color=colorEig[k])
		ax[fig].set_xlabel(nameData[choose[i][0]]);
		ax[fig].set_ylabel(nameData[choose[i][1]]);
		ax[fig].axis('equal')
		for item in ([ax[fig].title, ax[fig].xaxis.label, ax[fig].yaxis.label]):
			item.set_fontsize(fontSize);
		for item in (ax[fig].get_xticklabels() + ax[fig].get_yticklabels()):
			item.set_fontsize(fontSize-5);
		ax[fig].yaxis.set_ticks_position('left')
		ax[fig].xaxis.set_ticks_position('bottom')
		ax[fig].tick_params(length=10, width=5, which='major');
		ax[fig].tick_params(length=5, width=2, which='minor');

	fig+=1;
	# new data
	if len(pred)==0:
		ax[fig].plot(score[0,:],score[1,:],'o',color=colorPoint);
	else:
		col=["blue","red"];
		for i in xrange(len(pred)):
			ax[fig].plot(score[0,i],score[1,i],'o',color=col[pred[i]]);
	ax[fig].axis('equal')
	for item in ([ax[fig].title, ax[fig].xaxis.label, ax[fig].yaxis.label]):
		item.set_fontsize(fontSize);
	for item in (ax[fig].get_xticklabels() + ax[fig].get_yticklabels()):
		item.set_fontsize(fontSize-5);
	ax[fig].set_xlabel("$e_{1}$");
	ax[fig].set_ylabel("$e_{2}$");
	ax[fig].yaxis.set_ticks_position('left')
	ax[fig].xaxis.set_ticks_position('bottom')
	ax[fig].tick_params(length=10, width=5, which='major');
	ax[fig].tick_params(length=5, width=2, which='minor');

	#correlation

	# eigenvalue
	print(latent)
	fig+=1;
	ax[fig].bar([i for i in xrange(len(latent))],real(latent),width=0.5,color=colorEig);
	for item in ([ax[fig].title, ax[fig].xaxis.label, ax[fig].yaxis.label]):
		item.set_fontsize(fontSize);
	for item in (ax[fig].get_xticklabels() + ax[fig].get_yticklabels()):
		item.set_fontsize(fontSize-5);
	#ax[fig].set_xlabel("#eig");
	ax[fig].set_ylabel("$Re(\lambda)$");
	ax[fig].set_xticks(0.25+array([i for i in xrange(len(latent))]));
	ax[fig].set_xticklabels(["$\lambda_{%d}$"%(i+1) for i in xrange(len(latent))]);
	ax[fig].yaxis.set_ticks_position('left')
	ax[fig].xaxis.set_ticks_position('bottom')
	ax[fig].tick_params(length=10, width=5, which='major');
	ax[fig].tick_params(length=5, width=2, which='minor');

	#eigenvectors
	fig+=1;
	inter=2.;
	xMean=array([inter*k for k in xrange(len(coeff[0,:]))]);
	sizeBar=1./float(min(len(latent),nEigMax));
	for i in xrange(min(len(latent),nEigMax)):
		x=array([inter*k+(i/float(min(len(latent),nEigMax))) for k in xrange(len(coeff[0,:]))]);
		ax[fig].bar(x,coeff[i,:],width=0.2,color=colorEig[i]);
	ax[fig].plot([-0.5,inter*len(idx)],[0,0],'--',color="black",lw=3);
	ax[fig].set_xticks(xMean+0.5);
	ax[fig].set_xticklabels(nameData,rotation=90);
	ax[fig].set_xlim([-0.5,len(idx)*inter]);
	for item in ([ax[fig].title, ax[fig].xaxis.label, ax[fig].yaxis.label]):
		item.set_fontsize(fontSize);
	for item in (ax[fig].get_xticklabels() + ax[fig].get_yticklabels()):
		item.set_fontsize(fontSize-5);
	ax[fig].yaxis.set_ticks_position('left')
	ax[fig].xaxis.set_ticks_position('bottom')
	ax[fig].tick_params(length=10, width=5, which='major');
	ax[fig].tick_params(length=5, width=2, which='minor');
