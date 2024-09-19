from openpyxl import load_workbook
import numpy as np
import pandas as pd

class master_data:
    def __init__(self, dataPath):
        ### import data from excel workbook
        
        self.wb = load_workbook(dataPath)
        self.segWs = self.wb['SegmentProperties']
        self.cntWs = self.wb['BioProbeCountMatrix']
        
        self.segValues = [[y.value for y in x] for x in self.segWs[self.segWs.calculate_dimension()]]
        self.cntValues = [[y.value for y in x] for x in self.cntWs[self.cntWs.calculate_dimension()]]

        self.dropData = False
        
        
    def get_data(self):
        ### Convert nested list to a pandas dataFrame and extract expression data with labels
        cntData = self.cntValues
        cntCols = self.cntValues[0]
        df = pd.DataFrame(self.cntValues)
        cntIndex = [x[2] for x in self.cntValues[1:]]
        indexName = self.cntValues[0][2]
        cntDF = pd.DataFrame(self.cntValues[1:], index=cntIndex, columns=cntCols)
        cntDF.index.name = indexName
        self.counts = cntDF.iloc[:,12:]
        self.counts = self.counts.astype(np.float64)      # Convert datatype to float64
        self.priobeInfo = cntDF.iloc[:,:12]
        segCols = self.segValues[0]
        segIndex = [x[4].replace(' | ',('_')) for x in self.segValues[1:]]
        self.segData = pd.DataFrame(self.segValues[1:], index=segIndex, columns=segCols)

        self.dataOrig = self.counts.copy()
        self.dataLog1 = np.log2(self.counts+1)            # Log transform data for QC and analysis steps

        self.probeClass = False                           # Keep a copy of the original data before transformation or normalisation
        # self.probeClass = df.iloc[self.targIdx:,2]      ### Index needs updating here also
        # self.probeClass.rename(index=rowLabels, inplace=True)
        # self.probeClass.rename(index='ProbeClass', inplace=True)

        # ToDo: Need to update probeclass handling and labeling
        
        self.probeClassDict = {
            'Positive': 'A',
            'Negative': 'B',
            'Control': 'C',
            'Endogenous': 'E'}

        return(self.counts, self.priobeInfo, self.segData)

    def get_descriptors(self):
        pass

    def add_class_mean(self, df):
        pass

    def drop_AOIs(self, includes, writeOrig=False):
        pass
        
    def set_threshold(self, threshold):
        self.threshold = threshold
        # ToDo: Check that all values in master data are also included in threshold dataFrame
        # ToDo: Convert threshold data to 0/1 data if needed

    def drop_probes(self, labels):
        pass
        
