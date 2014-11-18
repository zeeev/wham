#!/usr/bin/python
import argparse, csv, os, sys, re #std python imports
import numpy as np
from sklearn.ensemble import RandomForestClassifier #RF classifier from SKlearn
from sklearn.cross_validation import cross_val_score #validation stats from SKlearn
import itertools
import multiprocessing as mp #allows for parallelization of the classification to speed up script.

#########################
#Args
#########################

parser=argparse.ArgumentParser(description="Runs RandomForest classifier on WHAM output VCF files to classify structural variant type. Appends WC and WP flags for user to explore structural variant calls. The output is a VCF file written to standard out.")
parser.add_argument("VCF", type=str, help="User supplied VCF with WHAM variants; VCF needs AT field data")
parser.add_argument("training_matrix", type=str, help="training dataset for classifier derived from simulated read dataset")
parser.add_argument("--filter", type=str, help="optional arg for filtering type one of : ['sensitive', 'specific']; defaults to output all data if filtering if argument is not supplied.")
parser.add_argument("--proc", type=str, help="optional arg for number of proceses to run with classifier; higher thread number will increase speed of classifier; defaults to 1")
arg=parser.parse_args()


#########################
#Functions
#########################

#class object for processing VCF files. 
class vcf:
	"""
	class vcf generates an iterator for looping through a vcf
	file. Can add various functionalities here to further process
	the vcf in a number of ways
	chunksize = number of lines to process at once for parallel applications. 
	"""
	def __init__(self,file):
		self.f = open(file,'r')
		#proces all of the header lines of the vcf. 
		header = True #boolean to continye looping through header
		info_boolean = False #need a boolean to trigger and append new INFO fields
		while header:
			self.line = self.f.readline()
			line = self.line.strip()
			line = line.split("\t") #split line on tabs
			if line[0][0] == '#': #process header lines
				if re.search("##FORMAT", line[0]) and info_boolean == False: #first instance of ##FORMAT..
					#need to append new INFO fields for the corresponding data
					print '##INFO=<ID=WC,Number=1,Type=String,Description="WHAM classifier vairant type">'
					print '##INFO=<ID=WP,Number=4,Type=Float,Description="WHAM probability estimate for each structural variant classification from RandomForest model">'
					info_boolean = True #reset boolean to 
				print "\t".join( line ) #print results to stdout 
			else:
				header = False #break out of the loop
		
	def __iter__(self):
		return self
	
	def next(self, chunksize=5000):
		cnt = 0 #boolean for chunking.
		return_array = [] #initialize empty array to store line data. 
		#check here if we are currently on last line, and raise StopIteration to exit next()
		if len(self.line) == 0: #
			raise StopIteration
		while cnt < chunksize:
			line = self.line
			if len( line ) == 0:
				return( return_array ) 
				break #break out of loop because we are at last line in file. 
			else:
				return_array.append( line ) 
			self.line = self.f.readline()
			cnt += 1
		return( return_array )






#parse the targets for ML. converts text list of classified data
#into a numerical dataset with links to the classified names
def parse_targets( target ):
	"""
	target = list of factors to be turned into numerical classifiers. 
	for machine learning classifiction. ie. converts INR, DEL, INV, 
	DUP into integer factors
	"""
	target = np.array(target) #convert to np array for ease in processing
	names = np.unique( target ) #unique names of SV types (factors)

	#now iterate through and classify to an integer for SKlearn
	cnt = 0
	target_numerical = np.zeros( target.shape[0] ) #generate empty dataset
	for name in names:
		idx = np.where( name == target )
		target_numerical[ idx ] = cnt
		cnt += 1
	
	#setup return data structure
	RV = {'names': names, 'target': target_numerical}
	#can use index of SV type in 'names' to get text based variant
	#call from 'target', where they've been converted to integers.
	return( RV )


#method to run observed data through the trained model to output
#a vcf friendly return of classified variant call and the prediction
#probabilities for each call
def classify_data( _x, clf, names ):
	"""
	_x = pass the col 8 from vcf
	clf = machine learning object
	names = string names, zero indexed of variant calls. 
	"""
	_x = np.array(_x)
	class_idx = int( clf.predict(_x) )#predict classifier. can link back to dataset['target_names']
	prediction = names[ class_idx ] #lookup text based name for classification

	class_probs = clf.predict_proba(_x)[0] #gives weights for your predictions 1:target_names
	#convert back to text comma separated list
	class_str = ",".join( [ str(i) for i in class_probs ] )

	#parse the two data fields into a string so they can be appended to the vcf file. 
	return_str = "WC=" + prediction + ";WP=" + class_str 
	return( return_str )


#A general parser that takes the data in VCF flag field and parses it into a 
#dictionary data structure. Can then obtain whatever data needed by using
# RV['key']; ie. RV['GT'] ...
def parse_vcf_data( vdat ):
	"""
	vdat = string; column 8 from VCF file with INFO fields. 
	"""
	#start by parsing the vcf data into a dictionary structure
	#will be keyed by ["XX="] =  data 
	dict = {}
	vdat = vdat.split(";")
	for v in vdat:
		try:
			v = v.split("=")
		except:
			print "not valid VCF file"
		dict[ v[0] ] = v[1] #setup dict struct

	#return the dictionary structure data
	return( dict )


#takes vcf field data and runs various filtering specs.
def run_filters( vdat, filtering = None ):
	""" 
	vdat - dictionary of INFO field from VCF line
	filtering - dictionary of fields to be filtered; defaults to None
	Currently implemented for sensitive and specific. Can modify the
	filters to return False anytime you want to not report results based
	on filteirng criterion from the INFO field. 
	"""
	pass_filt = True #will remain true until we do not satisfy some criterion
	
	if filtering == None:
		return( pass_filt ) #break out early 

	#sensitive is very perimssive
	elif filtering == "sensitive":
		if int( vdat['NC'] ) < 2:
			pass_filt = False
			return( pass_filt )
		if pass_filt:
			return( pass_filt )

	#specific mapping is more restrictive on the filtering.
	elif filtering == "specific":
		if vdat['ED'] == 'nan':
			pass_filt = False
			return( pass_filt )
		BE = vdat['BE'].split(',') 
		if int(BE[-1]) < 2:
			pass_filt = False
			return( pass_filt )
		if int( vdat['NC'] ) < 3:
			pass_filt = False
			return( pass_filt )
		if pass_filt:
			return( pass_filt )

	#elif filtering == "user_defined":
	#	....

	else:
		raise ValueError('Not a valid --filter argumuent\n please try running with --help arg for instructions')


#fuction will process line information and classify variant for a line in VCF file.
def process_vcf( info ):
	"""
	pass izip object of line object and other needed vars
	info[0] = list of vcf lines from VCF object iterator. 
	info[1] = clf object
	info[2] = dataset dictionary
	info[3] = filter arg supplied by user
	"""
	#parse the args to function
	line_list = info[0] #list of lines from VCF obj
	clf = info[1] #randomForest object
	dataset = info[2] #dataset with class names 
	filter = info[3] #filter arg supplied by user

	#iterate over lines in the chunked data
	return_list = []
	for line in line_list:
		line = line.strip().split("\t")
		vdat = parse_vcf_data( line[7] ) #parse all of vcf appended data
		filter_bool = run_filters( vdat, filtering=filter ) #boolean of whether line info passes filters
	
		if filter_bool:
			_x = vdat[ 'AT' ].split(",") #create list from data in 'AT' field 
			results = classify_data( _x, clf, dataset['target_names'] )
		
			line[7] = line[7] + ";" + results #append data to correct vcf column
			#print "\t".join( line ) #print results to stdout
			print_line = "\t".join( line )
			return_list.append( print_line )
		else:
			return_list.append( None )
	#return the full list		
	return( return_list )



#########################
#MAIN
#########################



###########
#import and assign training data
###########
#all sklearn data will be in 2D array [ nsamples X nfeatures]

sys.stderr.write("processing training file... \n" )
#iterate over training file. select out the numerical and classifier data
data  = []
target = []
with open(arg.training_matrix) as t:
        for line in csv.reader(t,delimiter='\t'):
		if line[0][0] == "#": #add in this statemnt to print error if user supplies files in wrong order.
			raise ValueError('not a valid WHAM training file. perhaps you supplied arguments in the wrong order? \n please try running with --help arg for instructions') 
        	target.append( line[-1] ) #always have targets [classified SV] as last column
        	d = [ float(i) for i in line[0:-1] ]
        	data.append( d )


#populate the training dataset in sciKitLearn friendly structure. 
dataset = {} #empty data
dataset[ 'data' ] = np.array( data ) #all training data into 2-D array

#turn our target list into integers and return target names
target_parse = parse_targets( target )
dataset[ 'target' ] = np.array( target_parse['target'] )
dataset[ 'target_names' ] = np.array( target_parse['names'] )



###########
#random forest classification
###########

#setup inital params
clf = RandomForestClassifier( n_estimators=250 )
#run RFC on dataset with target classifiers; runs the model fit
clf = clf.fit( dataset['data'], dataset['target'] )


######
#run some sanity checks here. 
######
training_stats = clf.feature_importances_ #array of variable importances for model.

#print training stats to user
train_list = [ str(i) for i in training_stats ] #convert to str for printing to user.
sys.stderr.write("\t Training weights for RandomForest classifier \n\t N = %d training variables\n" %( len(train_list) ) ) 
sys.stderr.write("\t %s\n" %( "\t".join( train_list ) ) ) 

#need cross validation here. uses sklearn.cross_validation
scores = cross_val_score( clf, dataset['data'], dataset['target'] )
avg_val = scores.mean() * 100  #average cross validation levels
sys.stderr.write("\t results from cross validation:\n\t %f%s \n" %( avg_val, '%' ) )


######
#prediction and output
######

sys.stderr.write("processing VCF file through classifier... \n" ) 

#load VCF file into class obj
vcf_file = vcf(arg.VCF)

#parse the number of processes to enact
if arg.proc == None:
	proc_num = 1
else:
	proc_num = int( arg.proc )

###
#setup multiprocessing for the classification of SVs
###
p = mp.Pool( processes = proc_num )
results = p.imap(process_vcf, itertools.izip( vcf_file, itertools.repeat(clf), itertools.repeat(dataset), itertools.repeat(arg.filter) ) )

#iterate over the results and feed to stdout
for r in results:
	for rv in r: #iterate over the list of returned results
		if rv != None: #only print results that pass filtering specs. 
			print rv #write output to std out


#final output to std err that the run has finished. 
sys.stderr.write("...classifier finished \n" )



