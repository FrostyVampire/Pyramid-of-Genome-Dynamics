import numpy as np
from scipy.stats import wasserstein_distance

def getBlockCounts(blocks):
    # Initialize array with zeros (size matches genome size)
    blockCounts = [0] * abs(blocks[0][4])

    for block in blocks:
        # The length is the 5th element (index 4)
        blockCounts[abs(block[4]) - 1] += 1

    return blockCounts

# Find synteny blocks between 2 genomes
def findSynteny(genome1, genome2):
    genome1size = len(genome1)
    genome2size = len(genome2)
    # block = []
    blocks = []

    # If either genome is empty, skip
    if genome1size == 0 or genome2size == 0:
        return getBlockCounts(blocks)

    # Because the genomes are circular, move the last gene in the first genome to position 0 until it is no longer the start of a block
    push = 0
    j = 0
    while j < genome2size:
        if push < min(genome1size, genome2size) and j > -1 and genome2[j] == genome1[-1] and genome2[
            (j + 1) % genome2size] == genome1[0]:
            if push == 0:
                startPos2 = (j + 1) % genome2size
            push += 1
            genome1.insert(0, genome1.pop())
            j -= 2
        j += 1

    # If the genomes are the same, return early
    if min(genome1size, genome2size) == push:
        blocks.append([push, 0, genome1size - 1, startPos2, (startPos2 - 1) % genome2size])
        return getBlockCounts(blocks)

    # Iterate through the genomes to find synteny blocks
    blocksTemp = []
    addedBlock = True
    while addedBlock:
        addedBlock = False
        i = 0
        while i < genome1size:
            length = 0
            longestBlockLength = 0

            # Skip deleted genes
            while i < genome1size and genome1[i] == -1:
                i += 1
            if i == genome1size:  # reached end of genome
                break

            # Compare with second genome (regular form)
            for j in range(genome2size):
                # First gene in block found, find block length
                if i < genome1size and genome1[i] != -1 and genome1[i] == genome2[j]:
                    length = 1
                    while (genome1[(i + length) % genome1size] != -1 and genome1[(i + length) % genome1size] == genome2[
                        (j + length) % genome2size]):
                        length += 1

                    # If the current block is longer than the previous longest block, replace it
                    if length > longestBlockLength:
                        inversed = False
                        longestBlockLength = length
                        startPos1 = i
                        startPos2 = j
                    j = j + length

            # Compare with second genome (inversed form)
            for j in range(genome2size - 1, 0, -1):
                # First gene in block found, find block length
                if i < genome1size and genome1[i] != -1 and genome1[i] == genome2[j]:
                    length = 1
                    while (length < genome1size and genome1[(i + length) % genome1size] != -1 and
                           genome1[(i + length) % genome1size] == genome2[(j - length) % genome2size]):
                        length += 1

                    # If the current block is longer than the previous longest block, replace it
                    if length > longestBlockLength:
                        inversed = True
                        longestBlockLength = length
                        startPos1 = i
                        startPos2 = j
                    j = j - length

            # Add the longest block found to the list of blocks
            if longestBlockLength > 0:
                addedBlock = True
                if inversed:
                    blocksTemp.append([(startPos1 - push) % genome1size, (startPos1 + longestBlockLength - push - 1),
                                       startPos2, (startPos2 - longestBlockLength + 1), -longestBlockLength])
                else:
                    blocksTemp.append([(startPos1 - push) % genome1size, (startPos1 + longestBlockLength - push - 1),
                                       startPos2, (startPos2 + longestBlockLength - 1), longestBlockLength])
                i += longestBlockLength
            # If no block found, delete the gene (because no block will ever be found for it)
            else:
                genome1[i] = -1
                i += 1

        if addedBlock:
            blocksTemp = sorted(blocksTemp, key=lambda x: abs(x[0]), reverse=True)
            clearOverlaps(blocksTemp)
            deleteGenes(genome1, genome2, blocksTemp)
            # Move all the new blocks to the list of all blocks
            while len(blocksTemp) != 0:
                blocks.append(blocksTemp.pop(0))

    blocks = sorted(blocks, key=lambda x: abs(x[4]), reverse=True)

    # Increase all indexes by 1 to start at 1 instead of 0
    # for i in range(len(blocks)):
    #     blocks[i][0] += 1
    #     blocks[i][1] += 1
    #     blocks[i][2] += 1
    #     blocks[i][3] += 1
    #     blocks[i][4] = abs(blocks[i][4])
    return blocks

# Find the frequency count of a given block count
def normalize(blockCount):
    # Ignore empty entries
    genomeSize = len(blockCount)
    relevantSize = 1 + next(i for i in range(genomeSize - 1, -1, -1) if blockCount[i] != 0)  # Find the first index from the back that is not 0

    # Calculate the total amount of blocks
    totalBlocks = 0
    for i in range(relevantSize):
        totalBlocks += blockCount[i]
    # print(totalBlocks)

    frequency = [0] * relevantSize
    if totalBlocks > 0:
        for i in range(relevantSize):
            frequency[i] = blockCount[i] / totalBlocks

    return frequency

def compareDistributions(distribution1, distribution2):
    # Find the shorter length between the two distributions
    length = min(len(distribution1), len(distribution2))

    # Coordinates must match the sliced length
    blockLengths = list(range(1, length + 1))

    u_weights = distribution1[:length]
    v_weights = distribution2[:length]
    sumU = sum(u_weights)
    sumV = sum(v_weights)
    if sumU <= 0 or sumV <= 0 or np.isnan(sumU) or np.isnan(sumV) or np.isinf(sumU) or np.isinf(sumV):
        return 99999

    return wasserstein_distance(
        u_values=blockLengths,
        v_values=blockLengths,
        u_weights=u_weights,
        v_weights=v_weights
    )


# Support functions for findSynteny
# Remove any overlapping blocks
def clearOverlaps(blocks):
    length = len(blocks)
    i = 0
    while i < length:
        j = i + 1
        while j < length:
            # Block j begins after block i and ends before block i, or it ends after block i begins and ends before i ends
            if doesOverlap(blocks[i][0], blocks[i][1], blocks[j][0], blocks[j][1]) or \
                    doesOverlap(blocks[i][2], blocks[i][3], blocks[j][2], blocks[j][3]):
                del blocks[j]
            else:
                j += 1
        i += 1

# Check if there is overlap between blocks i and j
def doesOverlap(starti, endi, startj, endj):
    return (starti <= max(startj, endj)) and (endi >= min(startj, endj))

# Delete all genes found in blocks from the genomes
def deleteGenes(genome1, genome2, blocks):
    for i in range(len(blocks)):
        for j in range(blocks[i][0], blocks[i][1] + 1):
            genome1[j % len(genome1)] = -1
        for j in range(min(blocks[i][2], blocks[i][3]), max(blocks[i][2], blocks[i][3]) + 1):
            genome2[j % len(genome2)] = -1