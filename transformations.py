#-----------------------------------------------------------------------------------
#@filename: transformations.py
#@author: Jeffrey Creighton
#@date:2/11/2013
#@date_modified: 10/31/2014
#@Description: This module contains the transformations used in Eureka
#   Notable changes: as each transform is needed, a check is added if( i != '') to
#   ensure that there aren't any spaces in the data. This might be problematic when
#   a transformation which relies on the length approaches.
#------------------------------------------------------------------------------------
#@packages:
from sympy import *
from firstStirling import s
from secondStirling import S
from boustrophedonNumbers import E
from divisorsDictionary import d
from moebiusDictionary import mu
from BinomialDictionary import C
#-------------------------------------------------------------------------------

#n = Symbol('n', integer=True)
#---------------------------transformationT1()----------------------------------------
#@Param: list of integers
#@Return: list of integers
#@Description: the return value IS the param value, the identity of the sequence
#--------------------------------------------------------------------------------------
def transformationT1(sequence):
    transformed_list = []
    for i in sequence:
    	i = int(i)
        transformed_list.append(i)
    return transformed_list
#---------------------------------------------------------------------------------------

#---------------------------transformationT2()--------------------------------------------
#@Param: A list of integers
#@Return: A list of integers
#@Description: the sequence will be transformed to its partial sums
#----------------------------------------------------------------------------------------
def transformationT2(sequence):
    sum_val = 0
    transformed_list = []
    for i in sequence:
    	i = int(i)
        sum_val = sum_val + i
        transformed_list.append(sum_val)
    return transformed_list
#--------------------------------------------------------------------------------------------

#---------------------------------transformationT3()-----------------------------------------
#@Param: list of integers
#@Return: list of integers
#@Description: the sequence will be transformed to it's partial sum of squares
#---------------------------------------------------------------------------------------------
def transformationT3(sequence):
    transformed_list = []
    sum_val = 0
    for i in sequence:
        i = int(i)
        sq = i*i
        sum_val = sum_val + sq
        transformed_list.append(sum_val)
    return transformed_list
#--------------------------------------------------------------------------------------

#-----------------------------transformationT4()----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: the return value will be the param value transformed as the inverse
#partial sums of its binomial coefficients
#---------------------------------------------------------------------------------------
def transformationT4(sequence):
    transformed_list = []
    length = len(sequence)
    if(length>50):
	length = 50
    i = 0
    counter = 0
    n = 0
    sum_val = 0
    neg_val =  1
    while(counter < length):
	i=0
	k=0
	neg_val = 1
	sum_val =0
	while(k<=n):
	    if(n==k or k==0):
		bin = 1
	    elif(k==1 or (k == (n-1))):
		bin = n
	    else:
		key = str(n)+","+str(k)
		bin = int(C[key])
	    transform_val = neg_val * int(sequence[i]) * bin
	    sum_val = sum_val + transform_val
	    k = k+1
	    i = i+1
	    neg_val = neg_val * -1
	transformed_list.append(sum_val)
	n = n + 1
	counter = counter + 1
    return transformed_list
#---------------------------------------------------------------------------------------

#------------------------------transformationT5()-----------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: return value will be the partial sums of self convolution of the param
#value.
#-------------------------------------------------------------------------------------
def transformationT5(sequence):
    transformed_list = []
    L = len(sequence)
    for n in range(L):
	k = 0
	sum_val = 0
	while(k<=n):
	    transform_val = int(sequence[k])*int(sequence[n-k])
	    sum_val = sum_val + transform_val
	    k = k + 1
	transformed_list.append(sum_val)
    return transformed_list
#-------------------------------------------------------------------------------------

#-----------------------transformationT6()--------------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: The return val is the Linear Weighted Partial Sums of the param value.
#-------------------------------------------------------------------------------------
def transformationT6(sequence):
    transformed_list = []
    length = len(sequence)
    sum_val = 0
    k = 1
    for i in range(length):
        i = int(i)
        trans_val = (k)*(int(sequence[i]))
        sum_val = sum_val + trans_val
        transformed_list.append(sum_val)
        k = k + 1
    return transformed_list
#------------------------------------------------------------------------------------

#---------------------------transformationT7()--------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: the return value will be the param value transformed as the partial sums
#of its binomial coeffiecient values
#----------------------------------------------------------------------------------
def transformationT7(sequence):
    transformed_list = []
    length = len(sequence)
    if(length > 50):
	length = 50
    i = 0
    counter = 0
    n = 0
    sum_val = 0
    while(counter < length):
	i = 0
        k = 0
	sum_val = 0
	while(k<=n):
	    if(n==k or k==0):
	        bin = 1
	    elif(k==1 or (k ==(n-1))):
		bin = n
	    else:
		key = str(n)+","+str(k) 
		bin = int(C[key])
	    transform_val = int(sequence[i]) * bin
	    sum_val = sum_val + transform_val
	    k = k + 1
	    i = i + 1
	transformed_list.append(sum_val)
	n = n + 1
	counter = counter + 1
    return transformed_list
#-------------------------------------------------------------------------------

#----------------------------transformationT8()----------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the sequence will be replaced by their product of two
#consecutive elements
#----------------------------------------------------------------------------------------
def transformationT8(sequence):
    transformed_list = []
    L = (len(sequence))-1
    n_position = 0
    for i in range(L):
	if(sequence[i] != ''):
        	n_value = int(sequence[n_position])
        	n_plus_value = int(sequence[n_position + 1])
        	product_value = (n_value) * (n_plus_value)
        	transformed_list.append(product_value)
        	n_position = n_position + 1
    return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT9()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Cassini transformation
#----------------------------------------------------------------------------------------
def transformationT9(sequence):
    transformed_list = []
    L = (len(sequence))-2
    n_position = 1
    for i in range(L):
	if(sequence[i] != ''):
		i = int(i)	
        	n_minus_value = int(sequence[n_position - 1])
        	n_plus_value = int(sequence[n_position + 1])
        	n_value = int(sequence[n_position])
        	outside = (n_minus_value)*(n_plus_value)
        	inside = (n_value)*(n_value)
        	cass_val = outside - inside
        	transformed_list.append(cass_val)
        	n_position = n_position + 1
    return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT10()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the First Stirling transformation
#----------------------------------------------------------------------------------------
def transformationT10(sequence):
	transformed_list = []
	L = len(sequence)
	if(L > 50):
		L = 50
	transformed_list.append(sequence[0])
	for n in range(L):
	    if(n != 0):
	    	k = 1
	    	sum_val = 0 
	    	while(k<=n):
		    key = str(n) + "," + str(k)
		    trans_val = s[key] * int(sequence[k])
		    sum_val = sum_val + trans_val
		    k = k + 1
	        transformed_list.append(sum_val)
	return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT11()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Second Stirling transformation
#----------------------------------------------------------------------------------------
def transformationT11(sequence):
	transformed_list = []
	L = len(sequence)
	if(L > 50):
		L = 50
	transformed_list.append(sequence[0])
	for n in range(L):
	    if(n != 0):
	        k = 1
	        sum_val = 0
	        while(k<=n):
		    key = str(n) + "," + str(k)
		    trans_val = S[key] * int(sequence[k])
		    sum_val = sum_val + trans_val
		    k = k + 1
	        transformed_list.append(sum_val)
	return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT12()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Boustrophedon transformation
#----------------------------------------------------------------------------------------
def transformationT12(sequence):
	transformed_list = []
	L = len(sequence)
	if(L > 50):
		L = 50
	for n in range(L):
	    k = 0
	    sum_val = 0
	    while(k <= n):
		if(k==n or k==0):
		    binomial_value = 1
		elif(k==1 or (k == (n-1))):
		    binomial_value = n
		else:
		    key = str(n)+","+str(k)
		    binomial_value = int(C[key])
		sequence_value = int(sequence[k])
		boust_value = E[((n)-k)]
		trans_value = binomial_value * boust_value * sequence_value
		sum_val = sum_val + trans_value
		k = k + 1
	    transformed_list.append(sum_val)
	return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT13()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the First Differences transformation
#----------------------------------------------------------------------------------------
def transformationT13(sequence):
	transformed_list = []
	n = 0
	for i in range(len(sequence)-1):
		first = int(sequence[n + 1])
		second = int(sequence[n])
		trans_val = first - second
		transformed_list.append(trans_val)
		n = n + 1
	return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT14()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Catalan transformation
#----------------------------------------------------------------------------------------
def transformationT14(sequence):
    transformed_list = []
    if(len(sequence) > 0 and len(sequence) < 400):
        transformed_list.append(int(sequence[0]))
        counter = 1
        n = 1
        length = len(sequence)-1
        #transformed_list.append(int(sequence[0]))
        while(n <= length):
	    sum_val = 0
	    k = 0 
	    while(k <= n):
	        bin = int(C[str(2*n-k-1)+","+str(n-k)])
	        trans_val = (k * bin * int(sequence[k]))
	        k = k + 1
    	        sum_val = sum_val + trans_val
	    sum_val = sum_val / n
	    n = n + 1
            transformed_list.append(sum_val)
    return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT15()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Hankel Transformation
#----------------------------------------------------------------------------------------
def transformationT15(sequence):
    transformed_list = []
    matrix_list = []

    transform_val = 0
    L = len(sequence)
    L = (L+1)/2

    for i in range(L+1):
        j = 0
        k = 0
        l = j + i
	
        m = Matrix([[sequence[(k+j)] for j in range(i)] for k in range(l)])
	if(i == 0):
	    i = i
	else:
	    transform_val = m.det()
	    if(len(str(transform_val)) > 100):
	        break
	    else:
	        transformed_list.append(transform_val)
	k = k + 1
	j = j + 1

    return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT16()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes  the Sum of Divisors transformation
#----------------------------------------------------------------------------------------
def transformationT16(sequence):
	transformed_list = []
	L = len(sequence)
	if(L > 50):
	    L = 50
	for n in range(L):
		sum_val = 0
		k = 1
		while(k <= n):
			sum_val = 0
			j = 0	
			while(j < k):
				key = str(k) + "," + str(j)
				b_value = d[key]
				a_value = int(sequence[b_value])
				sum_val = sum_val + a_value
				if(b_value == k):
					j = k
				else:
					j = j + 1
			k = k + 1
		transformed_list.append(sum_val)	
	return transformed_list
#------------------------------------------------------------------------------------------

#----------------------------transformationT17()-----------------------------------------
#@Param: a list of integers
#@Return: a list of integers
#@Description: each element of the list undergoes the Moebius transformation
#----------------------------------------------------------------------------------------
def transformationT17(sequence):
	transformed_list = []
	L = len(sequence)
	
	for n in range(L):
		k = 0
		sum_val = 0
		while(k <= n):
			sum_val = 0
			j = 0
			while(j <= k and k != 0):
				#might move  (k<0) into an if right here
				divisor_key = str(k) + "," + str(j)
				divisor_val = d[divisor_key]
				#if(sequence[divisor_val] != ''):
				s_val = int(sequence[divisor_val])
					
				mu_key = str(k/divisor_val)
				mu_val = mu[mu_key]
					
				transform_val = mu_val * s_val
				sum_val = sum_val + transform_val 
			
				if(divisor_val == k):
					j = k + 1
				else:
					j  = j + 1 
			k = k + 1
		transformed_list.append(sum_val)
	return transformed_list
#----------------------------------------------------------------------------------------
