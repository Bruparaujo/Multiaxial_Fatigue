# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 2017

@author: Bruno
"""


import math
import matplotlib.pyplot as plt
import numpy as np
import time
from MultiArray import MultiArray


def readFile(fileName):
    '''
    File must keep following order:
    [xStressAxis, yStressAxis, xyStressAxis, zStressAxis, xStrainAxis, yStrainAxis, xyStrainAxis, zStrainAxis]
    '''
    inputHistory = [[], [], [], [], [], [], [], []]
    numberOfColumns = len(inputHistory)
    with open(fileName, 'r') as f:
        for line in f:
            if line.startswith('stress'):
                continue
            else:
                splitStr = line.split(';')
                for column in range(numberOfColumns):
                    inputHistory[column].append(float(splitStr[column])) # lê cada coluna
    return inputHistory

def raceTrack(inputHistory, filterAmplitude, loadType):

    nDim = int(len(inputHistory) / 2)
    sqrtOf3 = 3**(1/2)

    if loadType == 0:
        getDim = 0
    elif loadType == 1:
        getDim = 4

    loadSize = len(inputHistory[0])     # "inputHistory[dim][line]"
    print("load size before racetrack = ", loadSize)

    filterAmplitudeSqr = filterAmplitude**2

    indexStack = []
    newIndexStack = [0, (loadSize-1)]

    '''partitioning step'''
    isNewKeyPoint = True
    while isNewKeyPoint == True:

        isNewKeyPoint = False
        indexStack = newIndexStack      # "update old index stack"
        newIndexStack = [0]             # "clear new index stack"

        numKeyPoint = len(indexStack)
        for j in range(1 , numKeyPoint):     # "sweeps all current key points"
            lastKP = indexStack[j-1]
            nextKP = indexStack[j]

            vectorKP = [0] * nDim
            normVectorKPSum = 0
            for dim in range(nDim):     # vector between two key points
                if dim == 2:
                    if loadType == 0:
                        lastLoad = inputHistory[dim + getDim][lastKP] * sqrtOf3
                        nextLoad  = inputHistory[dim + getDim][nextKP] * sqrtOf3
                    elif loadType == 1:
                        lastLoad = inputHistory[dim + getDim][lastKP] / sqrtOf3
                        nextLoad  = inputHistory[dim + getDim][nextKP] / sqrtOf3
                else:
                    lastLoad = inputHistory[dim + getDim][lastKP]
                    nextLoad  = inputHistory[dim + getDim][nextKP]

                vectorKP[dim] = nextLoad - lastLoad
                normVectorKPSum += ( vectorKP[dim] ) ** 2

            normVectorKP = normVectorKPSum ** (1/2)

            unitVectorKP = [0] * nDim
            if normVectorKP > 0:
                for dim in range(nDim):     # unit vector between two key points
                    unitVectorKP[dim] = vectorKP[dim] / normVectorKP
            else:
                print('NOTE: key points have same coordinate in partitioning step')
                pass

            distanceSqrMax = 0
            for index in range( (lastKP + 1) , nextKP ):      # sweep point between P[j-1]:P[j]

                sampleVector = [0] * nDim
                sampleProjection = 0
                sampleNormSum = 0
                for dim in range(nDim):     # vector between sample point and key point
                    if dim == 2:
                        if loadType == 0:
                            lastLoad = inputHistory[dim + getDim][lastKP] * sqrtOf3
                            sampleLoad = inputHistory[dim + getDim][index] * sqrtOf3
                        elif loadType == 1:
                            lastLoad = inputHistory[dim + getDim][lastKP] / sqrtOf3
                            sampleLoad  = inputHistory[dim + getDim][index] / sqrtOf3
                    else:
                        lastLoad = inputHistory[dim + getDim][lastKP]
                        sampleLoad = inputHistory[dim + getDim][index]

                    sampleVector[dim] = sampleLoad - lastLoad
                    sampleProjection += unitVectorKP[dim] * sampleVector[dim]    # compute the square of distance from sample point to key point segment
                    sampleNormSum += sampleVector[dim] ** 2    # compute the square of distance from sample point to key point segment

                distanceSqr = sampleNormSum - sampleProjection ** 2

                if distanceSqr > distanceSqrMax:    # keep track of max distant point
                    distanceSqrMax = distanceSqr
                    indexMax = index

            if distanceSqrMax > filterAmplitudeSqr:     # point is not filtered
                newIndexStack.append(indexMax)
                isNewKeyPoint = True

            newIndexStack.append(indexStack[j])     # add this key point to new index stack

    print("remaining points after partitioning step = ", len(newIndexStack))

    '''racetrack step'''
    indexStack = newIndexStack      # "update old index stack"
    newIndexStack = []             # "clear new index stack"
    numKeyPoint = len(indexStack)
    for j in range(1 , numKeyPoint):     # "sweeps all current key points"

        lastKP = indexStack[j-1]
        nextKP = indexStack[j]

        vectorKP = [0] * nDim
        normVectorKPSum = 0
        for dim in range(nDim):     # vector between two key points
            if dim == 2:
                if loadType == 0:
                    lastLoad = inputHistory[dim + getDim][lastKP] * sqrtOf3
                    nextLoad  = inputHistory[dim + getDim][nextKP] * sqrtOf3
                elif loadType == 1:
                    lastLoad = inputHistory[dim + getDim][lastKP] / sqrtOf3
                    nextLoad  = inputHistory[dim + getDim][nextKP] / sqrtOf3
            else:
                lastLoad = inputHistory[dim + getDim][lastKP]
                nextLoad  = inputHistory[dim + getDim][nextKP]

            vectorKP[dim] = nextLoad - lastLoad
            normVectorKPSum += ( vectorKP[dim] ) ** 2

        normVectorKP = normVectorKPSum ** (1/2)

        unitVectorKP = [0] * nDim
        if normVectorKP > 0:
            for dim in range(nDim):     # unit vector between two key points
                unitVectorKP[dim] = vectorKP[dim] / normVectorKP
        else:
            print('NOTE: key points have same coordinate in racetrack step')
            pass

        translationDirection = 0    # initial translation direction = 0
        lastIndex = lastKP
        sphereCenter = 0    # initial slotted plat center = 0
        for index in range( (lastKP + 1) , (nextKP + 1) ):      # sweep point between P[j-1]:P[j]+1

            sampleProjection = 0
            for dim in range(nDim):     # projection of sample point to key point segment
                if dim == 2: # make shear_xy to efective shear_xy
                    if loadType == 0:
                        lastLoad = inputHistory[dim + getDim][lastKP] * sqrtOf3
                        sampleLoad  = inputHistory[dim + getDim][index] * sqrtOf3
                    elif loadType == 1:
                        lastLoad = inputHistory[dim + getDim][lastKP] / sqrtOf3
                        sampleLoad  = inputHistory[dim + getDim][index] / sqrtOf3
                else:
                    lastLoad = inputHistory[dim + getDim][lastKP]
                    sampleLoad = inputHistory[dim + getDim][index]

                sampleProjection += ( sampleLoad - lastLoad ) * unitVectorKP[dim]

            distanceCenter = sampleProjection - sphereCenter
            distanceSqr = distanceCenter * distanceCenter

            if distanceSqr > filterAmplitudeSqr:

                reversedDirection = distanceCenter * translationDirection
                if reversedDirection <= 0:     # point is not filtered
                    newIndexStack.append(lastIndex)

                    if distanceCenter < 0:    # get sign of distanceCenter
                        translationDirection = -1
                    elif distanceCenter > 0:
                        translationDirection = 1
                    else:
                        translationDirection = 0

                sphereCenter +=  (abs(distanceCenter) - filterAmplitude) * translationDirection
                lastIndex = index

        if lastIndex < nextKP:
            newIndexStack.append(nextKP)
    newIndexStack.append(lastIndex)

    print("remaining points after racetrack step = ", len(newIndexStack))

    filteredData = []
    for dim in range(nDim*2):
        filteredData.append([])
    for index in newIndexStack:
        for dim in range(nDim*2):
            filteredData[dim].append(inputHistory[dim][index])

    return filteredData

def planeProject(inputHistory, modelToExecute, planeToExecute, tetaToExecute):

    xStress, yStress, xyStress, zStress, xStrain, yStrain, xyStrain, zStrain = inputHistory
    lenData = len(xStress)

    '''Teta convert degree to rad'''
    tetaInRad = (tetaToExecute * math.pi) / 180.
    cos2Teta = math.cos(2 * tetaInRad)
    cosTetaSqr = math.cos(tetaInRad) ** 2
    sin2Teta = math.sin(2 * tetaInRad)
    sinTetaSqr = math.sin(tetaInRad) ** 2

    shearStressA = []
    normalStress = []
    shearStressA45 = []
    shearStressB = []
    normalStress45 = []
    shearStrainA = []
    normalStrain = []
    shearStrainA45 = []
    shearStrainB = []
    normalStrain45 = []

    for line in range(lenData):
        '''Stress'''
        shearStressAPoint = xyStress[line] * cos2Teta + ( (yStress[line] - xStress[line]) * sin2Teta ) / 2.
        normalStressPoint = xyStress[line] * sin2Teta + xStress[line] * cosTetaSqr + yStress[line] * sinTetaSqr
        shearStressA45Point = shearStressAPoint / (2. ** 0.5)
        shearStressBPoint = (normalStressPoint + zStress[line]) / 2.
        normalStress45Point = (normalStressPoint - zStress[line]) / 2.

        '''Strain'''
        shearStrainAPoint = xyStrain[line] * cos2Teta + (yStrain[line] - xStrain[line]) * sin2Teta
        normalStrainPoint = xyStrain[line] * sin2Teta  / 2. + xStrain[line] * cosTetaSqr + yStrain[line] * sinTetaSqr
        shearStrainA45Point = shearStrainAPoint / (2. ** 0.5)
        shearStrainBPoint = normalStrainPoint - zStrain[line]
        normalStrain45Point = (normalStrainPoint + zStrain[line]) / 2.

        shearStressA.append(shearStressAPoint)
        normalStress.append(normalStressPoint)
        shearStressA45.append(shearStressA45Point)
        shearStressB.append(shearStressBPoint)
        normalStress45.append(normalStress45Point)
        shearStrainA.append(shearStrainAPoint)
        normalStrain.append(normalStrainPoint)
        shearStrainA45.append(shearStrainA45Point)
        shearStrainB.append(shearStrainBPoint)
        normalStrain45.append(normalStrain45Point)

    # Case A(T) plane
    if planeToExecute == 'AT':

        if modelToExecute == 'Goodman' or modelToExecute == 'ESWT':
            mainLoad = MultiArray(normalStress)
            auxiliarLoad = MultiArray(normalStress)

        elif modelToExecute == 'SWT':
            mainLoad = MultiArray(normalStrain)
            auxiliarLoad = MultiArray(normalStress)

        elif modelToExecute == 'LiuI':
            mainLoad = MultiArray(normalStrain)
            auxiliarLoad = MultiArray(normalStress, shearStrainA, shearStressA)

    # Case A(S) plane
    elif planeToExecute == 'A':

        if modelToExecute == 'Findley' or modelToExecute == 'McDiarmid':
            mainLoad = MultiArray(shearStressA)
            auxiliarLoad = MultiArray(normalStress)

        elif modelToExecute == 'Fatemi-Socie':
            mainLoad = MultiArray(shearStrainA)
            auxiliarLoad = MultiArray(normalStress)

        elif modelToExecute == 'Brown-Miller':
            mainLoad = MultiArray(shearStrainA)
            auxiliarLoad = MultiArray(normalStrain, normalStress)

        elif modelToExecute == 'LiuII':
            mainLoad = MultiArray(shearStrainA)
            auxiliarLoad = MultiArray(shearStressA, normalStrain, normalStress)

    # Case BII plane
    elif planeToExecute == 'B':

        if modelToExecute == 'Findley' or modelToExecute == 'McDiarmid':
            mainLoad = MultiArray(shearStressB)
            auxiliarLoad = MultiArray(normalStress45)

        elif modelToExecute == 'Fatemi-Socie':
            mainLoad = MultiArray(shearStrainB)
            auxiliarLoad = MultiArray(normalStress45)

        elif modelToExecute == 'Brown-Miller':
            mainLoad = MultiArray(shearStrainB)
            auxiliarLoad = MultiArray(normalStrain45, normalStress45)

        elif modelToExecute == 'LiuII':
            mainLoad = MultiArray(shearStrainB)
            auxiliarLoad = MultiArray(shearStressB, normalStrain45, normalStress45)

    # Case BII-III plane
    elif planeToExecute == 'BA':

        if modelToExecute == 'Findley' or modelToExecute == 'McDiarmid':
            mainLoad = MultiArray(shearStressB, shearStressA45)
            auxiliarLoad = MultiArray(normalStress45)

        elif modelToExecute == 'Fatemi-Socie':
            mainLoad = MultiArray(shearStrainB, shearStrainA45)
            auxiliarLoad = MultiArray(normalStress45)

        elif modelToExecute == 'Brown-Miller':
            mainLoad = MultiArray(shearStrainB, shearStrainA45)
            auxiliarLoad = MultiArray(normalStrain45, normalStress45)

        elif modelToExecute == 'LiuII':
            mainLoad = MultiArray(shearStrainB, shearStrainA45)
            shearStressMises = []
            for line in range(lenData):
                shearStressMises.append( (shearStressB[line]**2 + shearStressA45[line]**2 )**0.5 )
            auxiliarLoad = MultiArray(shearStressMises, normalStrain45, normalStress45)

    else:
        print('PLANO {} ! FORMATO INVÁLIDO'.format(planeToExecute))
        import sys
        sys.exit(-1)


    return mainLoad, auxiliarLoad

#####################################################################
'''Multiaxial rainflow - 1 load, N aux'''

def rainflow(mainLoads, auxLoads):
    'Sanity Check'
    assert mainLoads.count == 1, 'Esta funcao so aceita 1 main'
    assert isinstance(mainLoads, MultiArray), 'mainLoad deve ser uma Lista simples'
    assert isinstance(auxLoads, MultiArray), 'auxLoads deve ser uma Lista simples'

    # Inputs
    mainLoad = mainLoads.first()

    loadSize = len(mainLoad)

    # Temporary Stacks for active list
    mainStack = []
    # = MultiArray([], [], [])
    minAuxStacks = MultiArray(*( [] for _ in range(auxLoads.count) ))
    maxAuxStacks = MultiArray(*( [] for _ in range(auxLoads.count) ))


    # Outputs
    halfCycleRange = []
    # = MultiArray([], [], [])
    halfCycleValley = MultiArray(*( [] for _ in range(auxLoads.count) ))
    halfCyclePeak = MultiArray(*( [] for _ in range(auxLoads.count) ))


    for line in range(loadSize):

        # adiciona novo ponto ao final das listas ativas [mainStack, auxLoads, auxLoads]
        mainStack.append(mainLoad[line])
        for col in range(auxLoads.count):
            currentMinAuxStack = minAuxStacks.get(col)
            currentMaxAuxStack = maxAuxStacks.get(col)

            currentMinAuxStack.append(auxLoads.get(col)[line])
            currentMaxAuxStack.append(auxLoads.get(col)[line])

        # atualiza lista ativa: k -> k+=1
        k = len(mainStack)
        while k >= 3:

            if ((mainStack[-3] <= mainStack[-2]) and (mainStack[-2] <= mainStack[-1])) \
                    or ((mainStack[-3] >= mainStack[-2]) and (mainStack[-2] >= mainStack[-1])):     # No Peak or Valley or Stop point

                # avalia max e min Auxiliar
                for col in range(auxLoads.count):
                    currentMinAuxStack = minAuxStacks.get(col)
                    currentMaxAuxStack = maxAuxStacks.get(col)

                    currentMinAuxStack[-3] = min(currentMinAuxStack[-3], currentMinAuxStack[-2])
                    currentMaxAuxStack[-3] = max(currentMaxAuxStack[-3], currentMaxAuxStack[-2])

                # remove ponto não pico/vale (k-=1)
                del mainStack[-2]
                for col in range(auxLoads.count):
                    del minAuxStacks.get(col)[-2]
                    del maxAuxStacks.get(col)[-2]

            elif ((mainStack[-2] < mainStack[-3]) and (mainStack[-3] <= mainStack[-1])) \
                    or ((mainStack[-2] > mainStack[-3]) and (mainStack[-3] >= mainStack[-1])):    # critério (a) ou (b) - conta + 1/2 ciclo de S[-3] a S[-2]

                mainRange = abs(mainStack[-3] - mainStack[-2])
                halfCycleRange.append(mainRange)
                for col in range(auxLoads.count):
                    currentMinAuxStack = minAuxStacks.get(col)
                    currentMaxAuxStack = maxAuxStacks.get(col)

                    auxiliaryValley = min(currentMinAuxStack[-3], currentMinAuxStack[-2])
                    auxiliaryPeak = max(currentMaxAuxStack[-3], currentMaxAuxStack[-2])

                    outputMinAuxStack = halfCycleValley.get(col)
                    outputMaxAuxStack = halfCyclePeak.get(col)

                    outputMinAuxStack.append(auxiliaryValley)
                    outputMaxAuxStack.append(auxiliaryPeak)

                if k > 3:   # critério (b) - conta + 1/2 ciclo de S[-3] a S[-2]

                    halfCycleRange.append(mainRange)
                    for col in range(auxLoads.count):
                        currentMinAuxStack = minAuxStacks.get(col)
                        currentMaxAuxStack = maxAuxStacks.get(col)

                        auxiliaryValley = min(currentMinAuxStack[-3], currentMinAuxStack[-2])
                        auxiliaryPeak = max(currentMaxAuxStack[-3], currentMaxAuxStack[-2])

                        outputMinAuxStack = halfCycleValley.get(col)
                        outputMaxAuxStack = halfCyclePeak.get(col)

                        outputMinAuxStack.append(auxiliaryValley)
                        outputMaxAuxStack.append(auxiliaryPeak)

                    # retain max/min among last entries
                    for col in range(auxLoads.count):
                        currentMinAuxStack = minAuxStacks.get(col)
                        currentMaxAuxStack = maxAuxStacks.get(col)

                        currentMinAuxStack[-4] = min(currentMinAuxStack[-4], currentMinAuxStack[-3], currentMinAuxStack[-2])
                        currentMaxAuxStack[-4] = max(currentMaxAuxStack[-4], currentMaxAuxStack[-3], currentMaxAuxStack[-2])

                    # remove ciclo contado - k-=2
                    del mainStack[-3:-1]
                    for col in range(auxLoads.count):
                        del minAuxStacks.get(col)[-3:-1]
                        del maxAuxStacks.get(col)[-3:-1]

                else:       # critério (a) - remove 1/2 ciclo contado (k-=1)

                    del mainStack[-3]
                    for col in range(auxLoads.count):
                        del minAuxStacks.get(col)[-3]
                        del maxAuxStacks.get(col)[-3]

            else:
                break
            k = len(mainStack)

    # critério (c) - conta 1/2 ciclo de pares restantes
    k = len(mainStack)
    while k >= 2:
        # conta 1/2 ciclo de S[-1] a S[-2]
        mainRange = abs(mainStack[0] - mainStack[1])
        halfCycleRange.append(mainRange)
        for col in range(auxLoads.count):
            currentMinAuxStack = minAuxStacks.get(col)
            currentMaxAuxStack = maxAuxStacks.get(col)

            auxiliaryValley = min(currentMinAuxStack[0], currentMinAuxStack[1])
            auxiliaryPeak = max(currentMaxAuxStack[0], currentMaxAuxStack[1])

            outputMinAuxStack = halfCycleValley.get(col)
            outputMaxAuxStack = halfCyclePeak.get(col)

            outputMinAuxStack.append(auxiliaryValley)
            outputMaxAuxStack.append(auxiliaryPeak)

        # remove 1/2 ciclo contado e atualiza k (k-=1)
        del mainStack[0]
        for col in range(auxLoads.count):
            del minAuxStacks.get(col)[0]
            del maxAuxStacks.get(col)[0]
        k = len(mainStack)

    return halfCycleRange, halfCycleValley, halfCyclePeak

#####################################################################
'''Modified Wang-Brown - 2 loads, N aux'''

def normCalculator(thisEntry, nextEntry):
    '''
    Transformar em classe: calcula a norma entre 2 vetores
    '''

    thisEntryB = thisEntry[0]
    nextEntryB = nextEntry[0]
    thisEntryA = thisEntry[1]
    nextEntryA = nextEntry[1]

    normCalculated = math.sqrt( (nextEntryB - thisEntryB)**2 + (nextEntryA - thisEntryA)**2 )

    return normCalculated

def medianCalculator(thisEntry, nextEntry):

    medianCalculated = (nextEntry + thisEntry) / 2.

    return medianCalculated

def findMaxNorm(inputMWB):
    findMaxNorm_time = time.time()

    shearInput = inputMWB
    n = len(inputMWB)
    lastNorm = -1.

    lastI = 0
    lastJ = 1
    for i in range(1, n):
        for j in range(i+1, n):
            newNorm = normCalculator(shearInput[i], shearInput[j])
            if newNorm > lastNorm:
                lastNorm = newNorm
                lastI = i
                lastJ = j
            else:
                pass

    elapsed_time = time.time() - findMaxNorm_time
    print ('Time after for(i,j): {} s'.format(round(elapsed_time,3)))

    normI = normCalculator(shearInput[lastI], [0, 0])
    normJ = normCalculator(shearInput[lastJ], [0, 0])

    if normJ < normI:
        firstPoint = lastJ
    else:
        firstPoint = lastI

    outputMWB = inputMWB
    outputMWB.extend(inputMWB[0:firstPoint])
    del outputMWB[0:firstPoint]

    return outputMWB

def findAlfa(entryI, entryK1, entryJ, entryJ1):

    a = normCalculator(entryJ1, entryJ)
    if a == 0:
        print('Error finding α in findAlfa: division by zero')
        print(entryJ1, entryJ)
        exit()
    b = normCalculator(entryJ, entryI)
    c = normCalculator(entryJ1, entryI)
    p = normCalculator(entryK1, entryI)

    a2 = a**2
    b2 = b**2
    c2 = c**2
    p2 = p**2

    alfaPartI = (a2 + b2 - c2)
    alfaPartII = ( alfaPartI**2 - 4*a2 * (b2 - p2) ) ** 0.5

    alfaPlus = ( alfaPartI + alfaPartII ) / (2*a2)
    alfaMinus = ( alfaPartI - alfaPartII ) / (2*a2)

    lowerAlfa = min(alfaPlus, alfaMinus)
    higherAlfa = max(alfaPlus, alfaMinus)
    if (lowerAlfa >= 0) and (lowerAlfa <= 1.00001):
        validAlfa = lowerAlfa
    else:
        if (higherAlfa >= 0) and (higherAlfa <= 1.00001):
            validAlfa = higherAlfa
            if lowerAlfa > -0.00001:
                print('   ')
                print('lowerAlfa = {}'.format(lowerAlfa))
                print('   ')
        else:
            print('Error finding valid α: α- = {} α+ = {}'.format(alfaMinus, alfaPlus))
            exit()

    # if (alfaPlus >= 0) and (alfaPlus <= 1.00001):
    #     if (alfaMinus >= 0) and (alfaPlus <= 1.00001):
    #         if alfaMinus < alfaPlus:
    #             validAlfa = alfaMinus
    #         else:
    #             validAlfa = alfaPlus
    #     else:
    #         validAlfa = alfaPlus
    # elif (alfaMinus >= 0) and (alfaPlus <= 1.00001):
    #     validAlfa = alfaMinus
    # else:
    #     print('Error finding valid α: α- = {} α+ = {}'.format(alfaMinus, alfaPlus))
    #     exit()
    return validAlfa

def reorderCount(outMWB, lastpoint):

    lenOutMWB = len(outMWB)
    orderedCount=[]
    for i in range(lenOutMWB):
        firstOut = outMWB[i][0]

        isCloseLoop = False
        j=i
        while isCloseLoop == False and j < (lenOutMWB - 1):
            j+=1

            lastOut = outMWB[j][-1]
            if firstOut == lastOut:
                isCloseLoop = True
                print('''

                Loop inside MWB!

                ''')

        if isCloseLoop == False:
            if i < (lenOutMWB - 1):
                orderedCount.append(outMWB[i])
        else:
            firstHalf = MultiArray(outMWB[i])
            lastHalf = MultiArray(outMWB[j])

            if firstHalf.get(0)[-1] == lastHalf.get(0)[0]:
                del lastHalf.get(0)[0]
                fullLoop =  ( firstHalf.get(0) + lastHalf.get(0) )
                orderedCount.append(fullLoop)
            else:
                fullLoop =  ( firstHalf.get(0) + lastHalf.get(0) )
                orderedCount.append(fullLoop)

    return orderedCount

def MWB(mainLoads, auxLoads):
    'Sanity Check'
    assert mainLoads.count == 2, 'Esta funcao so aceita 2 main'
    assert isinstance(mainLoads, MultiArray), 'mainLoad deve ser uma Lista simples'
    assert isinstance(auxLoads, MultiArray), 'mainLoad deve ser uma Lista simples'

    mainLoad1 = mainLoads.get(0)
    mainLoad2 = mainLoads.get(1)

    historyInput = []
    for line in range( len(mainLoad1) ):
        if line == 0:
            auxStacks = MultiArray(*( [] for _ in range(auxLoads.count) ))
            for col in range(auxLoads.count):
                currentAuxStack = auxStacks.get(col)
                currentAuxStack.append(auxLoads.get(col)[line])

            historyInput.append( [mainLoad1[line], mainLoad2[line], auxStacks] )
        else:
            localAmpFilter = 0.001 /100.
            if abs(mainLoad1[line] - mainLoad1[line-1]) < localAmpFilter  and abs(mainLoad2[line] - mainLoad2[line-1]) < localAmpFilter:
                pass
            else:
                auxStacks = MultiArray(*( [] for _ in range(auxLoads.count) ))
                for col in range(auxLoads.count):
                    currentAuxStack = auxStacks.get(col)
                    currentAuxStack.append(auxLoads.get(col)[line])

                historyInput.append( [mainLoad1[line], mainLoad2[line], auxStacks] )

    loadInput = findMaxNorm(historyInput)

    shearB=0
    shearA=1
    auxLoad=2

    isClosedCycle = False
    if isClosedCycle == True:
        loadInput.append( loadInput[0] )          # Apenas para ciclo fechado
    else:
        pass

    lenLoadInput = len(loadInput) - 1
    alfa = [-1] * lenLoadInput
    count = []
    lastPoint = []
    lineI = 0
    while lineI <= (lenLoadInput-1):
        count.append( [loadInput[lineI]] )
        lastPoint.append(lineI)
        if alfa[lineI] > 0:
            tempCountB = loadInput[lineI][shearB] + alfa[lineI] * (loadInput[lineI+1][shearB] - loadInput[lineI][shearB])
            tempCountA = loadInput[lineI][shearA] + alfa[lineI] * (loadInput[lineI+1][shearA] - loadInput[lineI][shearA])
            count[lineI].append( [tempCountB, tempCountA, loadInput[lineI][auxLoad]] )
            lastPoint[lineI] = lineI + alfa[lineI]
            alfa[lineI] = 0

        elif alfa[lineI] == 0:
            del loadInput[lineI]
            del count[lineI]
            del alfa[lineI]
            lineI -= 1
            lenLoadInput -= 1

        else:
            count[lineI].append( loadInput[lineI+1] )
            alfa[lineI] = 0
            lineK = lineI
            lineJ = lineI + 1

            while lineJ <= (lenLoadInput-1):

                normJ1 = normCalculator(loadInput[lineJ+1], loadInput[lineI])
                normK1 = normCalculator(loadInput[lineK+1], loadInput[lineI])
                if normJ1 >= normK1:
                    newAlfa = findAlfa(loadInput[lineI], loadInput[lineK+1], loadInput[lineJ], loadInput[lineJ+1])
                    if alfa[lineJ] >= 0:
                        if alfa[lineJ] > newAlfa:
                            if alfa[lineJ] - newAlfa < 0.00001:     # Evita parada desnecessária (e-05)
                                pass
                            else:
                                tempCountB = loadInput[lineJ][shearB] + newAlfa * (loadInput[lineJ+1][shearB] - loadInput[lineJ][shearB])
                                tempCountA = loadInput[lineJ][shearA] + newAlfa * (loadInput[lineJ+1][shearA] - loadInput[lineJ][shearA])
                                count[lineI].append( [tempCountB, tempCountA, loadInput[lineJ][auxLoad] ] )

                                tempCountB = loadInput[lineJ][shearB] + alfa[lineJ] * (loadInput[lineJ+1][shearB] - loadInput[lineJ][shearB])
                                tempCountA = loadInput[lineJ][shearA] + alfa[lineJ] * (loadInput[lineJ+1][shearA] - loadInput[lineJ][shearA])
                                count[lineI].append( [tempCountB, tempCountA, loadInput[lineJ][auxLoad] ] )
                                lastPoint[lineI] = lineJ + alfa[lineJ]
                                alfa[lineJ] = newAlfa
                        else:
                            pass

                        break

                    else:
                        alfa[lineJ] = newAlfa

                        if alfa[lineJ] >= 0.99999:     # Evita parada desnecessária (e-05)          # NOTA: colocar antes na forma de elif
                            alfa[lineJ] = -1.0

                        else:
                            tempCountB = loadInput[lineJ][shearB] + alfa[lineJ] * (loadInput[lineJ+1][shearB] - loadInput[lineJ][shearB])
                            tempCountA = loadInput[lineJ][shearA] + alfa[lineJ] * (loadInput[lineJ+1][shearA] - loadInput[lineJ][shearA])
                            count[lineI].append( [tempCountB, tempCountA, loadInput[lineJ][auxLoad] ] )
                            count[lineI].append( loadInput[lineJ+1] )
                            lastPoint[lineI] = lineJ+1

                        lineK = lineJ
                else:
                    pass

                lineJ += 1
        lineI += 1

    orderedCount = reorderCount(count, lastPoint)

    return orderedCount

def MOI(inputMWB):
    '''
    outShearMWB[numCycles][numLoadCounted][typeLoad]
    '''

    typeShearB=0
    typeShearA=1
    typeAux=2

    '''Outputs'''
    halfCycleRange = []
    auxSize = inputMWB[0][0][typeAux].count
    halfCycleValley = MultiArray(*( [] for _ in range(auxSize) ))
    halfCyclePeak = MultiArray(*( [] for _ in range(auxSize) ))

    numCycles = len(inputMWB)
    for posCycle in range(numCycles):

        p = 0
        sumShearB = 0
        sumShearA = 0

        '''Temporary Stacks for active list'''
        minAuxStacks = MultiArray(*( [] for _ in range(auxSize) ))
        maxAuxStacks = MultiArray(*( [] for _ in range(auxSize) ))
        for col in range(auxSize):
            firstMinAuxStack = inputMWB[posCycle][0][typeAux].get(col)[0]
            firstMaxAuxStack = inputMWB[posCycle][0][typeAux].get(col)[0]

            currentMinAuxStack = minAuxStacks.get(col)
            currentMaxAuxStack = maxAuxStacks.get(col)

            currentMinAuxStack.append(firstMinAuxStack)
            currentMaxAuxStack.append(firstMaxAuxStack)

        numLoadCounted = int(len(inputMWB[posCycle]) /2 )
        for posCount in range(numLoadCounted):

            dShearBA = normCalculator( inputMWB[posCycle][posCount*2] , inputMWB[posCycle][posCount*2+1] )
            shearBMi = medianCalculator( inputMWB[posCycle][posCount*2][typeShearB] , inputMWB[posCycle][posCount*2+1][typeShearB] )
            shearAMi = medianCalculator( inputMWB[posCycle][posCount*2][typeShearA] , inputMWB[posCycle][posCount*2+1][typeShearA] )

            for col in range(auxSize):
                currentMinAuxStack = minAuxStacks.get(col)[0]
                currentMaxAuxStack = maxAuxStacks.get(col)[0]

                thisAux = inputMWB[posCycle][posCount*2][typeAux].get(col)[0]
                nextAux = inputMWB[posCycle][posCount*2+1][typeAux].get(col)[0]

                minStack = min([currentMinAuxStack, thisAux, nextAux])
                maxStack = max([currentMaxAuxStack, thisAux, nextAux])

                minAuxStacks.get(col)[0] = minStack
                maxAuxStacks.get(col)[0] = maxStack

            p += dShearBA
            sumShearB += shearBMi * dShearBA
            sumShearA += shearAMi * dShearBA

        if p == 0:
            print('Error finding perimeter in MOI: division by zero')
            print(inputMWB[posCycle], posCycle)
            import sys
            sys.exit(-1)

        shearBM = sumShearB / p
        shearAM = sumShearA / p

        sumIP=0
        for posCount in range(numLoadCounted):

            dShearBA = normCalculator( inputMWB[posCycle][posCount*2] , inputMWB[posCycle][posCount*2+1] )
            shearBMi = medianCalculator( inputMWB[posCycle][posCount*2][typeShearB] , inputMWB[posCycle][posCount*2+1][typeShearB] )
            shearAMi = medianCalculator( inputMWB[posCycle][posCount*2][typeShearA] , inputMWB[posCycle][posCount*2+1][typeShearA] )

            sumIP += ( ((dShearBA**2) / 12) + (shearBMi - shearBM)**2 + (shearAMi - shearAM)**2 ) * dShearBA

        IP = sumIP / p
        equivalentShear = (IP * 12) ** 0.5


        halfCycleRange.append(equivalentShear)

        for col in range(auxSize):
            currentMinAuxStack = minAuxStacks.get(col)[0]
            currentMaxAuxStack = maxAuxStacks.get(col)[0]

            outputMinAuxStack = halfCycleValley.get(col)
            outputMaxAuxStack = halfCyclePeak.get(col)

            outputMinAuxStack.append(currentMinAuxStack)
            outputMaxAuxStack.append(currentMaxAuxStack)

    return halfCycleRange, halfCycleValley, halfCyclePeak

#####################################################################

def newtonRaphson(inputValue):

    beta = inputValue[0]
    B = inputValue[1]
    gama = inputValue[2]
    C = inputValue[3]
    delta = inputValue[4]
    errorAdmissible = inputValue[5]

    icognita = min( math.log(delta / beta) / B , math.log(delta / gama) / C )

    errorGoal = 1

    i=1
    while errorGoal > 0:
        expB = beta * math.exp( B * icognita )
        expC = gama * math.exp( C * icognita )

        icognita -= ( expB + expC ) / ( B * expB + C * expC ) * math.log( ( expB + expC ) / delta )

        if ( expB + expC - delta ) == 0:
            print ('Newton Raphson has division by 0')
            print (beta, expB, B, gama, expC, C )
            print (delta, icognita, i, (expB + expC - delta), errorGoal )
            errorGoal = 1
        else:
            errorGoal = ( expB *  ( 1 + errorAdmissible ) ** B + expC *  ( 1 + errorAdmissible)  ** C  - delta ) / ( expB + expC - delta )
        #print (errorGoal)
        i += 1
        if i > 20:
            print ( 'Newton Raphson do not converge: iteration = {}'.format(i) )
            print ( delta, icognita, i, (expB + expC - delta), errorGoal )
            break
    #print ( delta, icognita )
    return icognita

def fatigueModels(countedCycles, parameters, modelToExecute):

    mainRange = countedCycles[0]
    auxiliaryValley = countedCycles[1]
    auxiliaryPeak = countedCycles[2]

    '''input parameters'''
    basicParameters = parameters[0]
    fatigueLifeParameters = parameters[1]
    modelParameter = parameters[2]

    elasticModulus = basicParameters[0]
    poisson = basicParameters[1]
    ultimateTensileStrength = basicParameters[2]
    yieldCyclicStrength = basicParameters[3]
    fatigueLimit = basicParameters[4]

    elasticExponent = fatigueLifeParameters[0]
    elasticCoefficient = fatigueLifeParameters[1]
    plasticExponent = fatigueLifeParameters[2]
    plasticCoefficient = fatigueLifeParameters[3]

    alfaStressScaleFactor = modelParameter[0]
    betaShearFatigueLimit = modelParameter[1]


    lenCount = len(mainRange)
    cummulatedDamage = 0
    damageInTime = []

    # Critério de Goodman
    if modelToExecute == 'Goodman':

        normalStressRange = mainRange
        normalStressMax = auxiliaryPeak.get(0)

        for i in range(lenCount):
            damageParameter = (normalStressRange[i] /2) / ( 1 - ((normalStressMax[i] - (normalStressRange[i] /2)) / ultimateTensileStrength ))
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            elif normalStressMax[i] > (ultimateTensileStrength - 5):
                '''
                By pass for large plastic deformation
                '''
                damageInTime.append(cummulatedDamage)
                pass
            else:
                eventLife = ( ( damageParameter / elasticCoefficient) ** ( 1/ elasticExponent) )
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de ESWT
    elif modelToExecute == 'ESWT':

        normalStressRange = mainRange
        normalStressMax = auxiliaryPeak.get(0)

        for i in range(lenCount):
            damageParameter = normalStressRange[i] * normalStressMax[i] / 2.
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = ( ( (damageParameter**0.5) / elasticCoefficient ) ** ( 1/ (elasticExponent) ) )
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de Findley
    elif modelToExecute == 'Findley' or modelToExecute == 'McDiarmid':

        shearA = mainRange
        normalStressMax = auxiliaryPeak.get(0)

        for i in range(lenCount):
            damageParameter = shearA[i] / 2. + alfaStressScaleFactor * normalStressMax[i]
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = (( damageParameter * (fatigueLimit / (betaShearFatigueLimit * elasticCoefficient)) ) ** ( 1/ elasticExponent)) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de SWT
    elif modelToExecute == 'SWT':

        normalStrainRange = mainRange
        normalStressMax = auxiliaryPeak.get(0)

        beta = elasticCoefficient ** 2 / elasticModulus
        exponentB = elasticExponent * 2
        gama = elasticCoefficient * plasticCoefficient
        exponentC = plasticExponent + ( elasticExponent )
        errorAdmissible = 1e-1

        for i in range(lenCount):
            damageParameter = abs(normalStrainRange[i] * normalStressMax[i] / 2.)
            inputValue = [beta, exponentB, gama, exponentC, damageParameter, errorAdmissible]
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = math.exp( newtonRaphson(inputValue) ) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de Brown-Miller
    elif modelToExecute == 'Brown-Miller':

        shearStrainRange = mainRange
        normalStrainPeak = auxiliaryPeak.get(0)
        normalStrainValley = auxiliaryValley.get(0)

        normalStressPeak = auxiliaryPeak.get(1)
        normalStressValley = auxiliaryValley.get(1)

        beta1 = (1 + poisson) + (1 - poisson) * alfaStressScaleFactor
        beta2 = (1.5 + 0.5 * alfaStressScaleFactor)

        # beta = beta1 * (elasticCoefficient) / elasticModulus
        exponentB = elasticExponent
        gama = beta2 * plasticCoefficient
        exponentC = plasticExponent
        errorAdmissible = 1e-3

        for i in range(lenCount):
            normalStressMedian = (normalStressPeak[i] + normalStressValley[i]) / 2.
            beta = beta1 * abs(elasticCoefficient - 2*normalStressMedian) / elasticModulus

            normalStrainRange = abs(normalStrainPeak[i] - normalStrainValley[i]) #/ 2
            damageParameter = abs(shearStrainRange[i] / 2.) + alfaStressScaleFactor * normalStrainRange
            inputValue = [beta, exponentB, gama, exponentC, damageParameter, errorAdmissible]
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = math.exp( newtonRaphson(inputValue) ) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de Fatemi-Socie
    elif modelToExecute == 'Fatemi-Socie':

        shearStrainRange = mainRange
        normalStressMax = auxiliaryPeak.get(0)

        beta = elasticCoefficient / elasticModulus
        exponentB = elasticExponent
        gama = plasticCoefficient
        exponentC = plasticExponent
        errorAdmissible = 1e-3

        for i in range(lenCount):
            damageParameter = abs(shearStrainRange[i]) / 2. * ( 1 + alfaStressScaleFactor * normalStressMax[i] / yieldCyclicStrength )
            inputValue = [beta, exponentB, gama, exponentC, damageParameter, errorAdmissible]
            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = math.exp( newtonRaphson(inputValue) ) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de LiuI
    elif modelToExecute == 'LiuI':

        normalStrainRange = mainRange
        normalStressPeak = auxiliaryPeak.get(0)
        normalStressValley = auxiliaryValley.get(0)
        shearStrainPeak = auxiliaryPeak.get(1)
        shearStrainValley = auxiliaryValley.get(1)
        shearStressPeak = auxiliaryPeak.get(2)
        shearStressValley = auxiliaryValley.get(2)

        beta = elasticCoefficient * plasticCoefficient
        exponentB = elasticExponent + plasticExponent
        gama = (elasticCoefficient ** 2) / elasticModulus

        exponentC = 2 * elasticExponent
        errorAdmissible = 1e-3

        for i in range(lenCount):
            normalStressRatio =   normalStressValley[i] / normalStressPeak[i]
            normalStressRange = abs(normalStressPeak[i] - normalStressValley[i]) #/ 2.
            shearStrainRange  = abs( shearStrainPeak[i] - shearStrainValley[i])  #/ 2.
            shearStressRange  = abs( shearStressPeak[i] - shearStressValley[i])  #/ 2.

            '''R > 0.7 -> crack is full open'''
            if normalStressRatio > 0.7:
                normalStressRatio = 0.7
            else:
                pass

            multiplierFactor = (2 / (1 - normalStressRatio)) * (1 / 4.)
            damageParameter = ( normalStrainRange[i] * normalStressRange + shearStrainRange * shearStressRange ) * multiplierFactor
            inputValue = [beta, exponentB, gama, exponentC, damageParameter, errorAdmissible]

            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = math.exp( newtonRaphson(inputValue) ) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    # Critério de LiuII
    elif modelToExecute == 'LiuII':

        shearStrainRange = mainRange
        shearStressPeak = auxiliaryPeak.get(0)
        shearStressValley = auxiliaryValley.get(0)
        normalStrainPeak = auxiliaryPeak.get(1)
        normalStrainValley = auxiliaryValley.get(1)
        normalStressPeak = auxiliaryPeak.get(2)
        normalStressValley = auxiliaryValley.get(2)

        beta = elasticCoefficient * plasticCoefficient
        exponentB = elasticExponent + plasticExponent
        gama = (elasticCoefficient ** 2) / elasticModulus
        exponentC = 2 * elasticExponent
        errorAdmissible = 1e-3

        for i in range(lenCount):
            normalStrainRange  = abs(normalStrainPeak[i] - normalStrainValley[i]) #/ 2.
            normalStressRange  = abs(normalStressPeak[i] - normalStressValley[i]) #/ 2.
            shearStressRange   = abs( shearStressPeak[i] -  shearStressValley[i]) #/ 2.

            normalStressMedian = (normalStressPeak[i] + normalStressValley[i]) / 2.
            normalElasticCoefficient = ( elasticCoefficient * math.sqrt(3) )
            divFactor = (normalElasticCoefficient - normalStressMedian)

            multiplierFactor = ( normalElasticCoefficient / divFactor )  * (1 / 4.)
            damageParameter = ( normalStrainRange * normalStressRange + shearStrainRange[i] * shearStressRange ) * multiplierFactor
            inputValue = [beta, exponentB, gama, exponentC, damageParameter, errorAdmissible]

            if damageParameter < 0:
                damageInTime.append(cummulatedDamage)
            else:
                eventLife = math.exp( newtonRaphson(inputValue) ) /2.
                eventDamage = 1. / ( eventLife * 2.0 )
                cummulatedDamage += eventDamage
                damageInTime.append(cummulatedDamage)

    else:
        print('Model to execute {} is not valid'.format(modelToExecute))
        import sys
        sys.exit(-1)

    return cummulatedDamage, damageInTime

def computationModulus(inputHistory, modelToExecute, planeToExecute, parameters, deltaPlaneToExecute):

    allPlanesLife=[]
    for candidatePlane in np.arange(0. ,180. , deltaPlaneToExecute):

        projectedMainLoad, projectedAuxLoad = planeProject(inputHistory, modelToExecute, planeToExecute, candidatePlane)

        if planeToExecute == 'BA':
            print(candidatePlane)
            outMWB = MWB(projectedMainLoad, projectedAuxLoad)
            countedLoad = MOI(outMWB)
        else:
            countedLoad = rainflow(projectedMainLoad, projectedAuxLoad)

        damageInPlane, damageInTime = fatigueModels(countedLoad, parameters, modelToExecute)
        allPlanesLife.append([damageInPlane, candidatePlane])

        if candidatePlane == 0:
            secondMinimumLife = [damageInPlane, candidatePlane]
            minimumLife = [damageInPlane, candidatePlane]
        elif minimumLife[0] < damageInPlane:
            secondMinimumLife[0] = minimumLife[0]
            secondMinimumLife[1] = minimumLife[1]
            minimumLife[0] = damageInPlane
            minimumLife[1] = candidatePlane
        else:
            if secondMinimumLife[0] < damageInPlane:
                secondMinimumLife[0] = damageInPlane
                secondMinimumLife[1] = candidatePlane
            else:
                pass
    print ("Accumulated damage for {} at plane {} is {} at {}° and {} at {}°".format(modelToExecute, planeToExecute, round(minimumLife[0],3), int(minimumLife[1]), round(secondMinimumLife[0],3), int(secondMinimumLife[1])))
    return allPlanesLife, damageInTime

def getPlotAxis(axis, column):
    plotYAxis = []
    axisSize = len(axis)
    for line in range(axisSize):
        plotYAxis.append(axis[line][column])
    return plotYAxis


if __name__ == "__main__":
    start_time = time.time()

    #####################################################################
    # Entries are: xAxis, yAxis, shearAxis or zAxis
    loadAxis = [ 'xAxis', 'yAxis', 'shearAxis', 'zAxis' ]
    # Entries are: STRESS or STRAIN
    loadingType = [ 'STRESS', 'STRAIN' ]
    # Entries are: AT, A, B or BA
    planeToExecute = [ 'AT', 'A', 'B', 'BA' ]
    # Entries are: Goodman, Findley or SWT
    modelToExecute = [ 'Goodman', 'Findley', 'SWT', 'Brown-Miller', 'ESWT', 'Fatemi-Socie', 'McDiarmid', 'LiuI', 'LiuII' ]
    #####################################################################

    '''Input history load'''
    filePath = 'E:\\User\Documentos\PUC\Fadiga\TCC\BKPs\BKP PENCIL 30.09.2019\PENCILDRIVE\MAFUI\src\services\diamond4.csv'
    fieldData = readFile(filePath)

    ##############################################################################################

    lineStyle = ['b', 'r', 'g', 'y']
    dotStyle = ['bo', 'r^', 'gs', 'yo']
    lowRange = 0
    highRange = 500

    '''Plotting input stress history'''
    for i in range(4):
        pointsToPlot = fieldData[i][lowRange:highRange]
        plt.plot(pointsToPlot, lineStyle[i], label=loadAxis[i])
        plt.plot(pointsToPlot, dotStyle[i])
    plt.title('Load')
    plt.xlabel('Points')
    plt.ylabel('Loads ' + loadingType[0])
    plt.legend()
    plt.show()

    plt.plot(fieldData[0][lowRange:highRange], fieldData[2][lowRange:highRange], lineStyle[0])
    plt.plot(fieldData[0][lowRange:highRange], fieldData[2][lowRange:highRange], dotStyle[0])
    plt.show()

    #####################################################################

    '''Plotting input strain history'''
    for i in range(4):
        pointsToPlot = fieldData[i+4][lowRange:highRange]
        plt.plot(pointsToPlot, lineStyle[i], label=loadAxis[i])
        plt.plot(pointsToPlot, dotStyle[i])
    plt.title('Load')
    plt.xlabel('Points')
    plt.ylabel('Loads ' + loadingType[0])
    plt.legend()
    plt.show()

    plt.plot(fieldData[4][lowRange:highRange], fieldData[6][lowRange:highRange], lineStyle[0])
    plt.plot(fieldData[4][lowRange:highRange], fieldData[6][lowRange:highRange], dotStyle[0])
    plt.show()

    #####################################################################

    '''Select to filter: stress or strain'''
    loadTtoFilter = 'strain'
    if loadTtoFilter == 'stress':
        loadType = 0
        filterAmplitude = 3     # 6%
        getDim = 0
    elif loadTtoFilter == 'strain':
        loadType = 1
        filterAmplitude = 0.00004     # 6%
        getDim = 4

    initInput = raceTrack(fieldData, filterAmplitude, loadType)

    ##############################################################################################

    lowRange = 0
    highRange = 500

    '''Plotting input stress history'''
    for i in range(4):
        pointsToPlot = initInput[i][lowRange:highRange]
        plt.plot(pointsToPlot, lineStyle[i], label=loadAxis[i])
        plt.plot(pointsToPlot, dotStyle[i])
    plt.title('Load')
    plt.xlabel('Points')
    plt.ylabel('Loads ' + loadingType[0])
    plt.legend()
    plt.show()

    plt.plot(initInput[0][lowRange:highRange], initInput[2][lowRange:highRange], lineStyle[1])
    plt.plot(initInput[0][lowRange:highRange], initInput[2][lowRange:highRange], dotStyle[0])
    plt.show()

    #####################################################################

    '''Plotting input strain history'''
    for i in range(4):
        pointsToPlot = initInput[i+4][lowRange:highRange]
        plt.plot(pointsToPlot, lineStyle[i], label=loadAxis[i])
        plt.plot(pointsToPlot, dotStyle[i])
    plt.title('Load')
    plt.xlabel('Points')
    plt.ylabel('Loads ' + loadingType[0])
    plt.legend()
    plt.show()

    plt.plot(initInput[4], initInput[6], lineStyle[1])
    plt.plot(initInput[4], initInput[6], dotStyle[0])
    plt.show()

    #####################################################################

    '''Parameters input'''

    '''material property'''
    elasticModulus = 193.0 * 10**3   # MPa                                          # E
    poisson = 0.3                                                                   # v
    shearElasticModulus = elasticModulus / (2*(1 + poisson))   # MPa                # G

    ultimateTensileStrength = 587.   # MPa                                          # S_U
    yieldCyclicStrength = 407   # MPa                                               # S_Yc

    if ultimateTensileStrength > 1400.:
        uniaxialFatigueLimit = 700.
    else:
        uniaxialFatigueLimit = ultimateTensileStrength / 2.   # MPa                 # S_L

    fatigueLimitRatio = 0.577
    torsionalFatigueLimit = uniaxialFatigueLimit * fatigueLimitRatio   # MPa        # TAU_L

    '''basquin property'''
    elasticExponent = -0.092                                                        # b
    elasticCoefficient = 745.4   # MPa                                              # SIGMA_c
    plasticExponent = -0.419                                                        # c
    plasticCoefficient = 0.161   # mm/mm                                            # EPSILON_c

    shearElasticExponent = elasticExponent                                          # b_gama
    shearElasticCoefficient = elasticCoefficient / math.sqrt(3)                     # TAU_c
    shearPlasticExponent = plasticExponent                                          # c_gama
    shearPlasticCoefficient = plasticCoefficient * math.sqrt(3)                     # GAMA_c

    '''default model property'''
    alfaStressScaleFactor = 0.5                                                     # default alfa
    betaShearFatigueLimit = torsionalFatigueLimit   # MPa                           # default beta

    '''default model property'''
    basicParameters = [elasticModulus, poisson, ultimateTensileStrength, yieldCyclicStrength, uniaxialFatigueLimit]
    shearBasicParameters = [shearElasticModulus, poisson, ultimateTensileStrength, yieldCyclicStrength, torsionalFatigueLimit]
    fatigueLifeParameters = [elasticExponent, elasticCoefficient, plasticExponent, plasticCoefficient]
    shearFatigueLifeParameters = [shearElasticExponent, shearElasticCoefficient, shearPlasticExponent, shearPlasticCoefficient]
    modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]

    #####################################################################
    computational_time = time.time()

    deltaPlane = 1   # candidate plane search precision

    # '''Damage evaluation for Goodman in plane A(T)'''
    # parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[0], planeToExecute[0], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[0]+' in plane '+planeToExecute[0])

    # '''Damage evaluation for Findley in plane A'''
    # alfaStressScaleFactor = (1 - uniaxialFatigueLimit/ ( 2 * torsionalFatigueLimit ) ) / ((uniaxialFatigueLimit/torsionalFatigueLimit - 1)**0.5)
    # betaShearFatigueLimit = (uniaxialFatigueLimit / 2.) / ((uniaxialFatigueLimit/torsionalFatigueLimit - 1)**0.5)
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[1], planeToExecute[1], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[1]+' in plane '+planeToExecute[1])

    # '''Damage evaluation for Findley in plane B'''
    # pulsatingFatigueLimit = uniaxialFatigueLimit * ultimateTensileStrength / (uniaxialFatigueLimit + ultimateTensileStrength)
    # alfaStressScaleFactor = (uniaxialFatigueLimit - pulsatingFatigueLimit) / (2*pulsatingFatigueLimit - uniaxialFatigueLimit)
    # betaShearFatigueLimit = (uniaxialFatigueLimit * pulsatingFatigueLimit / 2.) / (2*pulsatingFatigueLimit - uniaxialFatigueLimit)
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[1], planeToExecute[2], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[1]+' in plane '+planeToExecute[2])

    # '''Damage evaluation for McDiarmid in plane A'''
    # alfaStressScaleFactor = torsionalFatigueLimit / (2 * ultimateTensileStrength)
    # betaShearFatigueLimit = torsionalFatigueLimit
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[6], planeToExecute[1], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[6]+' in plane '+planeToExecute[1])

    # '''Damage evaluation for McDiarmid in plane B'''
    # torsionalFatigueLimitB = ( uniaxialFatigueLimit / 2 ) / ( 1 - ( uniaxialFatigueLimit / (4 * ultimateTensileStrength)))
    # alfaStressScaleFactor = torsionalFatigueLimitB / (2 *ultimateTensileStrength)
    # betaShearFatigueLimit = torsionalFatigueLimit
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[6], planeToExecute[2], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[6]+' in plane '+planeToExecute[2])

    # '''Damage evaluation for ESWT in plane A(T)'''
    # parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[4], planeToExecute[0], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[4]+' in plane '+planeToExecute[0])

    '''Damage evaluation for SWT in plane A(T)'''
    parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[2], planeToExecute[0], parameters, deltaPlane)
    plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[2]+' in plane '+planeToExecute[0])

    '''Damage evaluation for Brown-Miller in plane A'''
    alfaStressScaleFactor = 0.4                                         # alfa
    modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[3], planeToExecute[1], parameters, deltaPlane)
    plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[3]+' in plane '+planeToExecute[1])

    # '''Damage evaluation for Brown-Miller in plane B'''
    # alfaStressScaleFactor = 0.4                                         # alfa
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[3], planeToExecute[2], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[3]+' in plane '+planeToExecute[2])

    '''Damage evaluation for Fatemi-Socie in plane A'''
    alfaStressScaleFactor = 0.4     # alfa
    modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[5], planeToExecute[1], parameters, deltaPlane)
    plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[5]+' in plane '+planeToExecute[1])

    # '''Damage evaluation for Fatemi-Socie in plane B'''
    # alfaStressScaleFactor = 0.4     # alfa
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[5], planeToExecute[2], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[5]+' in plane '+planeToExecute[2])

    '''Damage evaluation for LiuI in plane A(T)'''
    parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[7], planeToExecute[0], parameters, deltaPlane)
    plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[7]+' in plane '+planeToExecute[0])

    # '''Damage evaluation for LiuII in plane A'''
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[8], planeToExecute[1], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[8]+' in plane '+planeToExecute[1])

    # '''Damage evaluation for LiuII in plane B'''
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[8], planeToExecute[2], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[8]+' in plane '+planeToExecute[2])

    #####################################################################
    '''
    Plane BA
    '''

    # '''Damage evaluation for Findley in plane BA'''
    # pulsatingFatigueLimit = uniaxialFatigueLimit * ultimateTensileStrength / (uniaxialFatigueLimit + ultimateTensileStrength)
    # alfaStressScaleFactor = (uniaxialFatigueLimit - pulsatingFatigueLimit) / (2*pulsatingFatigueLimit - uniaxialFatigueLimit)
    # betaShearFatigueLimit = (uniaxialFatigueLimit * pulsatingFatigueLimit / 2.) / (2*pulsatingFatigueLimit - uniaxialFatigueLimit)
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[1], planeToExecute[3], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), label=modelToExecute[1]+' in plane '+planeToExecute[3])

    # '''Damage evaluation for Brown-Miller in plane BA'''
    # alfaStressScaleFactor = 0.4     # alfa
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [basicParameters, fatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[3], planeToExecute[3], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), dotStyle[1], label=modelToExecute[3]+' in plane '+planeToExecute[3])

    # '''Damage evaluation for Fatemi-Socie in plane BA'''
    # alfaStressScaleFactor = 0.4     # alfa
    # modelParameter = [alfaStressScaleFactor, betaShearFatigueLimit]
    # parameters = [shearBasicParameters, shearFatigueLifeParameters, modelParameter]
    # allPlanesLife, damageInTime = computationModulus(initInput, modelToExecute[5], planeToExecute[3], parameters, deltaPlane)
    # plt.plot(getPlotAxis(allPlanesLife, 0), dotStyle[0], label=modelToExecute[5]+' in plane '+planeToExecute[3])

    # elapsed_time = time.time() - computational_time
    # print ('Elapsed time is: {} s'.format(round(elapsed_time,3)))
    #####################################################################

    '''Plotting Accumulated damage per Candidate Plane'''
    plt.xlabel('Candidate Plane')
    plt.ylabel('Accumulated damage')
    plt.title('Damage per plane ')
    plt.grid(True)
    plt.legend()
    plt.show()

    #####################################################################

    elapsed_time = time.time() - start_time
    print ('Elapsed time is: {} s'.format(round(elapsed_time,3)))

