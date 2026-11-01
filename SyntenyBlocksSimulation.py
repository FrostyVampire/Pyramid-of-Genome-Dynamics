# The user inputs a genome length and edge length
# The program then simulates the reversal process with the input parameters and compares the results to the results give by a formula

import random
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import os
import SBLD

numpy.set_printoptions(legacy='1.25')

# Find median size of blocks
def findMedian(total: int, blockCount):
    genomesize = int(simGenomeSize.get())

    if total == 0:
        return 0

    median = total / 2.0
    # If the total is odd
    if total % 2 == 1:
        for i in range(genomesize):
            median -= blockCount[i]
            if median <= 0:
                return i + 1

    # Else the total is even
    for i in range(genomesize):
        median -= blockCount[i]
        if median < 0.500000000001:
            i += 1
            if median > -0.499999999999:
                for j in range(genomesize):
                    if blockCount[j] != 0:
                        return (j + 1 + i) / 2.0
            return i
    return -1

# Display results in GUI
def displayResults(blockCount):
    genomeSize = int(simGenomeSize.get())
    relevantSize = 1 + next(i for i in range(genomeSize - 1, -1, -1) if blockCount[i] != 0)  # Find the first index from the back that is not 0

    # Simulation parameters from input
    jumpRate = float(edgeLengthJump.get()) if jumpCheckVar.get() else 0
    reversalRate = float(edgeLengthReversal.get()) if reversalCheckVar.get() else 0
    insRate = float(insertionRate.get()) if indelCheckVar.get() else 0.0
    delRate = float(deletionRate.get()) if indelCheckVar.get() else 0.0
    jumpSegLength = int(segmentLengthJump.get()) if jumpCheckVar.get() else 0
    indelSegLen = int(segmentLengthIndel.get()) if indelCheckVar.get() else 0

    # Calculate the total amount of blocks
    totalBlocks = 0
    for i in range(genomeSize):
        totalBlocks += blockCount[i]
    # print(totalBlocks)

    # Find block frequency
    frequency = SBLD.normalize(blockCount)

    # Find average
    avg = 0.0
    if totalBlocks > 0:
        for i in range(genomeSize):
            avg += blockCount[i] * (i + 1)
        avg /= totalBlocks

    for widget in frameRight.winfo_children():
        widget.destroy()

    # Labels for average and median
    label_avg = tk.Label(frameRight, text=f"Average Block Length: {avg:.2f}", font=("Arial", 12))
    label_avg.pack(pady=5)
    label_median = tk.Label(frameRight, text=f"Median Block Length: {findMedian(totalBlocks, blockCount)}",
                            font=("Arial", 12))
    label_median.pack(pady=5)

    # Initiate arrays and Q/P values
    formulaDistribution = [0] * genomeSize  # Probability that a maximal block of length k starts at a specific position
    formulaBlocks = [0] * relevantSize  # Estimated block count using the formula

    Qjump = numpy.exp(-jumpRate)  # No jump
    Pjump = 1 - Qjump             # Jump occurs

    Qrev = numpy.exp(-reversalRate)  # No reversal
    Prev = 1 - Qrev                  # Reversal

    phi = 1/(1+((insRate-delRate)/2))
    Qins = 1 / (1 + insRate)          # No insertion at slot
    Qdel = numpy.exp(-phi * delRate)  # No deletion
    Pins = 1 - Qins                   # Insertion occurs
    Pdel = 1 - Qdel                   # Deletion occurs
    P = Pins + Pdel - (Pins * Pdel)


    if jumpCheckVar.get():
        # k > seglen
        for k in range(jumpSegLength, relevantSize+1):
            formulaDistribution[k - 1] = 4*(Qjump**(2*k+jumpSegLength-2))*(Pjump*Pjump)
        # k == seglen
        formulaDistribution[jumpSegLength-1] += Pjump*(Qjump**(3*jumpSegLength-3))
        # k < seglen
        for k in range(jumpSegLength - 1, 0, -1):
            formulaDistribution[k - 1] = (1-Qjump**3)**2 * (Qjump**(3*k-3))

    # Formula (Reversal)
    elif reversalCheckVar.get():
        for k in range(1, relevantSize+1):
            formulaDistribution[k-1] = 4 * (Qrev ** (2*k-1)) * (Prev ** 2)

    # Formula (Indel)
    elif indelCheckVar.get():
        for k in range(1, relevantSize+1):
            formulaDistribution[k-1] = (Qins ** (k-1)) * (Qdel ** k) * (P ** 2)

    # Count blocks
    for k in range(1, relevantSize + 1):
        formulaBlocks[k-1] = formulaDistribution[k-1] * genomeSize

    # Normalize total blocks
    times = int(simTimes.get())
    for k in range(relevantSize+1):
        blockCount[k] = blockCount[k] / times

    print(f"Simulation parameters: {jumpRate} ({jumpSegLength}) | {reversalRate} | {insRate}, {delRate} ({indelSegLen}) | size {genomeSize}")
    trimmedBlocks = numpy.trim_zeros(blockCount, 'b')
    print(trimmedBlocks)
    # print("Formula distribution:")
    # print(formulaDistribution)

    # Find Wasserstein distance between (frequency, formulaDistribution)
    label_wasserstein = tk.Label(frameRight,
                                 text=f"Wasserstein Distance: {SBLD.compareDistributions(frequency, SBLD.normalize(formulaDistribution))}",
                                 font=("Arial", 12))
    label_wasserstein.pack(pady=5)

    # Prepare data for histogram
    lengths = [(i + 1) for i in range(genomeSize) if blockCount[i] > 0]
    valuesSim = [blockCount[i] for i in range(relevantSize+1) if blockCount[i] > 0]
    valuesFormula = [formulaBlocks[i] for i in range(relevantSize+1) if blockCount[i] > 0]

    # Plot with matplotlib
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lengths, valuesFormula, color='red', label='Formula', marker='s')  # Red line with squares
    ax.plot(lengths, valuesSim, color='blue', label='Simulation', marker='o')  # Blue line with circles
    ax.set_title("Synteny block distribution")
    ax.set_xlabel("Synteny Block Length")
    ax.set_ylabel("Block Count")

    # Save results. 1 file for graph and 1 file for values
    outputDir = "outputs"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    basefilename = f"jump{str(jumpRate).replace('.', '-')}_rev{str(reversalRate).replace('.', '-')}_ins{str(insRate).replace('.', '-')}_del{str(delRate).replace('.', '-')}_{genomeSize}_x{times}"
    # Check if file already exists
    filename = basefilename
    counter = 1
    while os.path.exists(os.path.join(outputDir, f"{filename}.png")):
        filename = f"{basefilename}({counter})"
        counter += 1
    ax.legend()
    graphPath = os.path.join(outputDir, f"{filename}.png")
    fig.savefig(graphPath)

    txtPath = os.path.join(outputDir, f"{filename}.txt")
    with open(txtPath, "w", encoding="utf-8") as f:
        f.write(f"Parameters: Size={genomeSize}, Jump={jumpRate} ({jumpSegLength}) | Reversal={reversalRate} | Insertion={insRate}, Deletion={delRate} ({indelSegLen}) | Times={times}\n")
        f.write(f"Qjump: {Qjump:.6f}, Qrev: {Qrev:.6f}, Qins: {Qins:.6f}, Qdel: {Qdel:.6f}\n")
        f.write(f"{'=' * 60}\n")
        f.write(f"{'Length':<8}\t{'Simulated':<12}\t{'Formula':<12}\t{'Ratio':<10}\n")
        f.write(f"{'-' * 60}\n")

        for i in range(len(trimmedBlocks)):
            simVal = blockCount[i]
            formVal = formulaBlocks[i] if i < len(formulaBlocks) else 0
            ratio = simVal / formVal if formVal > 0 else 0
            f.write(f"{i + 1:<8}\t{simVal:<12.4f}\t{formVal:<12.4f}\t{ratio:<10.4f}\n")

    # Embed plot into Tkinter
    canvas = FigureCanvasTkAgg(fig, master=frameRight)
    canvas.draw()
    canvas.get_tk_widget().pack()

    # Print results - How many blocks of each size and how far off are we.
    print("Printing results")
    for i in range(jumpSegLength+5):
        if blockCount[i] > 0:
            print(str(i+1) + ": " + str(blockCount[i]) + " / " + str(formulaBlocks[i]) + " = " + str(blockCount[i] / formulaBlocks[i]))
    print()

# Print all blocks
def printBlocks(blocks):
    for block in blocks:
        print(str(block[0]) + ": " + str(block[1]) + "-" + str(block[1] + block[3] - 1) + ", " + str(
            block[2]) + "-" + str(block[2] + block[3] - 1))

# Reversal model - take genome and edge length and perform reversals to create a new genome
def mutateGenome(genome, jumpRate, reversalRate, insRate, delRate, jumpSegLen, indelSegLen):
    #print("\n")
    #print(genome)
    genomesize = len(genome)

    # Determine how many events should occur based on a poisson distribution
    jumpCount = numpy.random.poisson(jumpRate * genomesize)
    reversalCount = numpy.random.poisson(reversalRate * genomesize)
    insCount = numpy.random.poisson(insRate * genomesize)
    delCount = numpy.random.poisson(delRate * genomesize)
    totalEvents = jumpCount + reversalCount + insCount + delCount
    eventNames = ["Jump", "Reversal", "Insertion", "Deletion"]

    # Run the mutations to create a new genome
    output = genome.copy()
    for i in range(totalEvents):
        # Choose a gene
        gene1 = random.randint(0, genomesize - 1)

        # Randomize between jump, reversal, insertion or deletion
        match random.choices(eventNames, weights=[jumpCount, reversalCount, insCount, delCount], k=1)[0]:
            case "Jump":
                #seglen = numpy.random.poisson(jumpSegLen) # We can safely assume that seglen << genomesize
                seglen = jumpSegLen
                if seglen == 0:
                    continue
                # Get segment
                if gene1 + seglen <= genomesize:
                    segment = output[gene1 : gene1 + seglen]
                    del output[gene1: gene1 + seglen]
                else:
                    headLen = seglen - (genomesize - gene1)
                    segment = output[gene1:] + output[:headLen]
                    del output[gene1:]
                    del output[:headLen]
                #print(f"Segment {segment} jumps")

                targetPos = random.randint(0, len(output))           # Choose the jump position
                output[targetPos:targetPos] = segment                   # Translocate the segment to a random spot
                jumpCount -= 1

            case "Reversal":
                # Choose the second gene
                gene2 = random.randint(0, genomesize - 1)
                if gene1 > gene2:   # Make sure gene1 is to the left of gene2
                    gene1, gene2 = gene2, gene1

                j = 0
                while j <= (gene2-gene1)/2:
                    output[gene1+j], output[gene2-j] = output[gene2-j], output[gene1+j]
                    j += 1
                reversalCount -= 1

            case "Insertion":
                seglen = max(1, numpy.random.poisson(indelSegLen))
                for i in range(seglen):
                    output.insert(gene1, -1) # Insert a new gene "-1" in the slot
                    genomesize += 1
                insCount -= 1

            case "Deletion":
                seglen = max(1, min(numpy.random.poisson(indelSegLen), len(output) - gene1))
                for i in range(seglen):
                    output.pop(gene1)   # Delete gene
                    genomesize -= 1
                delCount -= 1
    #print(output)
    return output

'''# Compare 2 distributions
def compareDistributions(distributionReal, distributionSimulation):
    # Sort the real distribution into bins as evenly as possible (max 10 bins)
    bins = [0]
    bestBins = bins
    binMaxLength = [0]  # The maximum block length (minus 1 because arrays start at 0) found in the bin i
    bestBinMaxLength = binMaxLength
    # Try 10 bins, then 9 bins, then 8... until 2 bins and find the number that gives the best ratio, give some leeway to bigger bin counts
    for binCount in range(10, 1, -1):
        # Reset the current array of bins
        bins = [0] * binCount
        binMaxLength = [-1] * binCount
        # Handle bins
        for i in range(0, binCount):
            for j in range(max(binMaxLength) + 1, len(distributionReal)):
                bins[i] += distributionReal[j]
                binMaxLength[i] = j
                if bins[i] >= 1 / binCount:
                    break

        if min(bestBins) == 0 or (min(bins) > 0 and max(bins) / min(bins) < (max(bestBins) / min(bestBins)) - 0.01):
            bestBins = bins
            bestBinMaxLength = binMaxLength

    # Remove bins of capacity 0
    for i in range(len(bestBins) - 1, 0, -1):
        if bestBins[i] == 0:
            bestBins.pop(i)
            bestBinMaxLength.pop(i)

    # Sort the simulated distribution into bins following the ranges we got from the real data
    binsSimulation = [0] * len(bestBins)
    for j in range(0, bestBinMaxLength[0] + 1):
        binsSimulation[0] += distributionSimulation[j]
    for i in range(1, len(bestBins)):
        for j in range(bestBinMaxLength[i - 1] + 1, bestBinMaxLength[i] + 1):
            binsSimulation[i] += distributionSimulation[j]
        # Add all remaining blocks to the last bin
        if i == len(bestBins) - 1:
            for j in range(bestBinMaxLength[i] + 1, len(distributionSimulation)):
                binsSimulation[i] += distributionSimulation[j]

    binCoordinates = list(range(len(bestBins)))

    return wasserstein_distance(
        u_values=binCoordinates,
        v_values=binCoordinates,
        u_weights=bestBins,
        v_weights=binsSimulation
    )
'''

# Run the simulation many times and take the average
def runSimulation():
    genomesize = int(simGenomeSize.get())
    times = int(simTimes.get())
    totalBlockCount = [0] * genomesize

    # Simulation parameters from input
    jumpRate = float(edgeLengthJump.get()) if jumpCheckVar.get() else 0
    reversalRate = float(edgeLengthReversal.get()) if reversalCheckVar.get() else 0
    insRate = float(insertionRate.get()) if indelCheckVar.get() else 0.0
    delRate = float(deletionRate.get()) if indelCheckVar.get() else 0.0
    jumpSegLen = int(segmentLengthJump.get()) if jumpCheckVar.get() else 0
    indelSegLen = int(segmentLengthIndel.get()) if indelCheckVar.get() else 0


    for i in range(times):
        # Initiate the first genome
        genome1 = [0] * genomesize
        for j in range(genomesize):
            genome1[j] = j+1

        # Get the second genome using the input parameters
        genome2 = mutateGenome(genome1, jumpRate, reversalRate, insRate, delRate, jumpSegLen, indelSegLen)

        # Compare the 2 genomes
        blockCount = SBLD.getBlockCounts(SBLD.findSynteny(genome1, genome2))
        for j in range(len(blockCount)):
            totalBlockCount[j] += blockCount[j]
    displayResults(totalBlockCount)

# GUI setup
# Toggle the state of entry boxes based on checkbox value.
def toggleInput(CheckVar, Entries):
    state = "normal" if CheckVar.get() else "disabled"
    for entry in Entries:
        entry.config(state=state)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Genome Reversal Synteny Simulation")
    root.geometry("900x700")

    frameLeft = tk.Frame(root, padx=10, pady=10)
    frameLeft.pack(side=tk.LEFT, fill=tk.Y)

    titleLabel = tk.Label(frameLeft, text="Genome Reversal Synteny Simulator", font=("Arial", 16))
    titleLabel.pack(pady=10)

    tableFrame = tk.Frame(frameLeft)
    tableFrame.pack(fill=tk.X, pady=10)

    tk.Label(tableFrame, text="Simulation Type", font=("Arial", 10, "bold")).grid(row=0, column=0, sticky="w", padx=5, pady=5)
    tk.Label(tableFrame, text="Lambda", font=("Arial", 10, "bold")).grid(row=0, column=1, padx=5, pady=5)
    tk.Label(tableFrame, text="Seg Length", font=("Arial", 10, "bold")).grid(row=0, column=2, padx=5, pady=5)

    # Jump
    jumpCheckVar = tk.BooleanVar(value=False)
    jumpCheck = tk.Checkbutton(tableFrame, text="Jump", variable=jumpCheckVar,
                               command=lambda: toggleInput(jumpCheckVar, [edgeLengthJump, segmentLengthJump]))
    jumpCheck.grid(row=1, column=0, sticky="w", padx=5, pady=5)

    edgeLengthJump = ttk.Entry(tableFrame, width=10)
    edgeLengthJump.insert(0, "0.2")
    edgeLengthJump.config(state="disabled")
    edgeLengthJump.grid(row=1, column=1, padx=5, pady=5)

    segmentLengthJump = ttk.Entry(tableFrame, width=10)
    segmentLengthJump.insert(0, "1")
    segmentLengthJump.config(state="disabled")
    segmentLengthJump.grid(row=1, column=2, padx=5, pady=5)


    # Reversal
    reversalCheckVar = tk.BooleanVar(value=False)
    revCheck = tk.Checkbutton(tableFrame, text="Reversal", variable=reversalCheckVar,
                              command=lambda: toggleInput(reversalCheckVar,[edgeLengthReversal]))
    # revCheck = tk.Checkbutton(tableFrame, text="Reversal", variable=reversalCheckVar,
    #                           command=lambda: toggleInput(reversalCheckVar,
    #                                                       [edgeLengthReversal, segmentLengthReversal]))
    revCheck.grid(row=2, column=0, sticky="w", padx=5, pady=5)

    edgeLengthReversal = ttk.Entry(tableFrame, width=10)
    edgeLengthReversal.insert(0, "0.2")
    edgeLengthReversal.config(state="disabled")
    edgeLengthReversal.grid(row=2, column=1, padx=5, pady=5)

    # segmentLengthReversal = ttk.Entry(tableFrame, width=10)
    # segmentLengthReversal.insert(0, "1")
    # segmentLengthReversal.config(state="disabled")
    # segmentLengthReversal.grid(row=2, column=2, padx=5, pady=5)

    # Indels
    indelCheckVar = tk.BooleanVar(value=False)
    indelCheck = tk.Checkbutton(tableFrame, text="Indels", variable=indelCheckVar,
                                command=lambda: toggleInput(indelCheckVar,
                                                            [insertionRate, deletionRate, segmentLengthIndel]))
    indelCheck.grid(row=3, column=0, sticky="w", padx=5, pady=5)

    # Inner container for Ins/Del alignment
    indelRatesFrame = tk.Frame(tableFrame)
    indelRatesFrame.grid(row=3, column=1, padx=5, pady=5)

    ttk.Label(indelRatesFrame, text="Ins:").pack(side=tk.LEFT)
    insertionRate = ttk.Entry(indelRatesFrame, width=6)
    insertionRate.insert(0, "0.05")
    insertionRate.config(state="disabled")
    insertionRate.pack(side=tk.LEFT, padx=(0, 5))

    ttk.Label(indelRatesFrame, text="Del:").pack(side=tk.LEFT)
    deletionRate = ttk.Entry(indelRatesFrame, width=6)
    deletionRate.insert(0, "0.05")
    deletionRate.config(state="disabled")
    deletionRate.pack(side=tk.LEFT)

    segmentLengthIndel = ttk.Entry(tableFrame, width=10)
    segmentLengthIndel.insert(0, "1")
    segmentLengthIndel.config(state="disabled")
    segmentLengthIndel.grid(row=3, column=2, padx=5, pady=5)

    # Global
    ttk.Separator(frameLeft, orient='horizontal').pack(fill='x', pady=20)

    ttk.Label(frameLeft, text="Number of simulations").pack()
    simTimes = ttk.Entry(frameLeft)
    simTimes.insert(0, "1")
    simTimes.pack()

    ttk.Label(frameLeft, text="Genome size (simulation)").pack()
    simGenomeSize = ttk.Entry(frameLeft)
    simGenomeSize.insert(0, "5000")
    simGenomeSize.pack()

    runButton = ttk.Button(frameLeft, text="Run Simulation", command=runSimulation)
    runButton.pack(pady=20)

    # Note
    color_guide = tk.Label(frameLeft, text="Red = formula\nBlue = induced")
    color_guide.pack(pady=10)

    frameRight = tk.Frame(root, padx=10, pady=10)
    frameRight.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    root.mainloop()