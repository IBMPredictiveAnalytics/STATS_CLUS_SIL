
#/***********************************************************************
# * Licensed Materials - Property of IBM
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2020
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp.
# ************************************************************************/

# Compute and plot silhouettes for clustering

__author__ = "SPSS, JKP"
__version__ = "1.0.1"

# history
# 06-sep-2011  original version
# 01-dec-2022  handle case where cluster assignment is missing

helptext="""
Calculate Cluster silhouette statistics for a clustering analysis already carried out.
This command requires the data in memory and takes time that is proportional
to the square of the number of cases.  With larger datasets, you may want to carry
this analysis out on a random sample of the data.

STATS CLUS SIL CLUSTER = varname VARIABLES = list of clustering variables
DISSIMILARITY = {EUCLID*|ABSDIFF|MINKOWSKI|GOWER|MAXIMUM}
[MINKOWSKIPOWER = number]
NEXTBEST = varname SILHOUETTE = varname
[/OPTIONS [RENUMBERORDINAL={NO*|YES}]
[VARWEIGHTS=numbers]
[MISSING={RESCALE*|IGNORE|OMIT}]]
[/OUTPUT [HISTOGRAM={YES*|NO}] [ORIENTATION={HORIZONTAL*|VERTICAL}]
[THREEDBAR={YES*|NO}] [THREEDCOUNTS={YES|NO*}]

Example:
STATS CLUS SIL CLUSTER=clusternumber VARIABLES=x1 x2 x3
DISSIMILARITY=ABSDIFF NEXTBEST = altclus SILHOUETTE = sil.

CLUSTER specifies the variable containing the cluster number for each case.
Gaps in the numbering, such as might occur if the cluster results are sampled,
are treated as clusters of size 0.

VARIABLES lists the variables to be used for calculating silhouettes.  These would
normally be the variables used to compute the clustering, but this is not required.

DISSIMILARITY specifies the function used to calculate the distance between cases.
EUCLID is standard Euclidean distance.
ABSDIFF is the sum of the absolute values of the variable differences.
MINKOWSKI is the pth root of the sum of the pth power of the absolute differences.
GOWER distinguishes variables according to the measurement level.  Scale variables
are measured by absolute difference with each difference divided by the variable
range.  Ordinal variables are renumbered so that the
values range from 1 to M, where M is the number of distinct values and then treated
as scale.  Nominal variable differences are 0 or 1 according to whether the values are
equal or not.
MAXIMUM is the maximum variable difference.  It is never rescaled, but the OMIT
setting is honored.

NEXTBEST and SILHOUETTE specify variables for the casewise output.
The NEXTBEST variable holds the number of the next best cluster.
The SILHOUETTE variable holds the silhouette statistics.

RENUMBERORDINAL causes ordinal variables to be recoded as above.  It is always
used for GOWER, but it can be used with any measure.

VARWEIGHTS can specify weights for each of the clustering variables.  If used,
each variable difference is multiplied by the corresponding weight.  As many weights
must be specified as there are clustering variables.

MISSING specifies how missing values are handled in case comparisions.  If a variable's
contribution to the measure cannot be computed because either case has a missing value,
there are three choices.  
RESCALE causes the measure to be scaled up in proportion to the number 
of contributions that are missing.  
IGNORE causes the measure to be used without any adjustment (as long as 
at least one variable contributes).  
OMIT causes the case to be omitted from the calculations if any variable's contribution
cannot be computed.

HISTOGRAM specifies whether a histogram of silhouette scores for each cluster
is produced.  By default, the histograms are panelled horizontally.  
ORIENTATION = VERTICAL causes the histograms to be stacked.  Since all the
histograms appear on one row, VERTICAL may be needed if there are many clusters.

THREEDBAR = YES produces a 3-d bar chart showing the average silhouette value by
assigned cluster and the next best clusters.

THREEDCOUNTS produces a 3-d bar chart showing for each assigned cluster
the percentage distribution of next best clusters.

/HELP displays this help and does nothing else.
"""

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd

import random, math

# debugging
        # makes debug apply only to the current thread
#try:
    #import wingdbstub
    #import threading
    #wingdbstub.Ensure()
    #wingdbstub.debugger.SetDebugThreads({threading.get_ident(): 1})
#except:
    #pass

class Dissimfuncs(object):
    """Dissimilarity calculators"""
    
    RESCALE, IGNORE, OMIT = [1,2,3]
    
    def __init__(self, cases, variables, dissim, missing, renumberordinal, weights, minkowskipower, vardict):
        """dissim is the dissimilarity measure to use
        vardict is a variable dictionary
        cases is the data matrix with first column containing the cluster number
        variables is the list of clustering variables
        dissim is the keyword name of the dissimilarity function to use
        if missing is rescale, dissimilarity measures are rescaled to compensate for variable
        differences that could not be computed because of missing values.  if missing is ignore, the 
        differences for those variables are just ignored.  If omit, no dissimilarity is calculated for
        that pair of cases.
        if renumberordinal, ordinal variable values are renumbered to be consecutive positive integers
        if weights, dissimilarities are weighted by the corresponding variable weight"""
        
        dissimfuncs={"euclid": self.euclid, "maximum" : self.maxdiff, "absdiff" : self.absdiff, "gower": self.gower,
            "minkowski" : self.minkowski}
        dissimlabels = {"euclid": _("Euclid"), "maximum" : _("Maximum"), "absdiff" : _("Absolute Difference"), 
            "gower": _("Gower"), "minkowski": "Minkowski(%s)"}
        
        self.cases = cases
        self.varmeta = [(vardict[v].VariableType, vardict[v].VariableLevel) for v in variables]
        if dissim in ["euclid", "maximum", "absdiff", "minkowski"]:  # these measures are numeric only
            if max([t for t,l in self.varmeta]) > 0:
                raise ValueError(_("""All variables must be numeric for the selected dissimilarity measure.""")) 
        self.dissimfunc = dissimfuncs[dissim]  # selected dissimilarity function
        self.dissimlabel = dissimlabels[dissim]
        if dissim == "minkowski":
            if minkowskipower is None:
                raise ValueError(_("The Minkowski measure requires the power to be specified"""))
            self.power = minkowskipower
            self.dissimlabel = dissimlabels[dissim] % minkowskipower
        self.intnumvars = len(cases[0])
        self.numvars = float(self.intnumvars)
        self.numclusvars = self.numvars - 1
        self.intnumclusvars =  int(self.numclusvars)
        self.numcases = len(self.cases)
        if missing == "rescale":
            self.missing = Dissimfuncs.RESCALE  # increase dissim in proportion to number of missing terms
        elif missing == "ignore":
            self.missing = Dissimfuncs.IGNORE  # ignore missing terms
        elif missing == "omit":
            self.missing = Dissimfuncs.OMIT  # no dissimilarity for cases with missing terms
            
        if dissim == "gower":
            renumberordinal = True
            
        if renumberordinal:
            self.renumber()
            
        # extra data prep required for gower measure
        # require ranges of predictors and rescaling of ordinal variables

        if dissim == "gower":
            # calculate ranges for clustering variables.  Ignore cluster number in first column
            mins = [min(item[k] for item in self.cases if item[k] is not None) for k in range(1, self.intnumvars)]
            maxes = [max(item[k] for item in self.cases if item[k] is not None) for k in range(1, self.intnumvars)]
            self.vrange = []
            try:
                for v in range(self.intnumclusvars):
                    if self.varmeta[v][1] != "nominal":   #numeric.  Nonnumeric variables can only be nominal
                        therange = maxes[v] - mins[v]
                        self.vrange.append(therange)
                    else:
                        self.vrange.append(None)
            except:
                raise ValueError(_("""At least one variable is entirely missing"""))
                     
        if weights:
            if len(weights) != self.numclusvars:
                raise ValueError(_("""The number of variable weights must be equal to the number of variables"""))
        self.weights = weights

    def renumber(self):
        """Renumber each ordinal variable so that its values are consective whole numbers starting at 1"""
        
        # The first column of self.cases is the cluster number and is excluded from the renumbering
        # varmeta excludes the cluster variable.
        
        for v in range(1, self.intnumvars):
            if self.varmeta[v - 1][1] == "ordinal":
                # construct sorted list of distinct,  nonmissing values
                thevalues = sorted(set([item[v] for item in self.cases]))
                if thevalues[0] is None:
                    del(thevalues[0])
                # recode ordinal values to 1:V where V is the number of distinct values
                recodedict = dict((k,float(v)) for k, v in zip(thevalues, list(range(1, len(thevalues)+1))))
                recodedict[None] = None
                for c in range(self.numcases):
                    self.cases[c][v] = recodedict[self.cases[c][v]]
        
    def calc(self, case1, case2):
        """calculate dissimilarity using specific function
        
        case1 and case2 are the observations whose distance will be calculated
        zeroth column is the cluster number and should be ignored"""
        
        # index identifies the variable in the clustering list
        # screen out variable values that can't be compared
        z = [(item, index) for index, item in enumerate(zip(case1[1:], case2[1:])) if item[0] is not None and item[1] is not None]
        lenz = len(z)
        
        # if only comparing complete cases or all values are missing, return None
        if (self.missing == Dissimfuncs.OMIT and lenz != self.intnumclusvars) or lenz == 0:
            return None
        
        #z consists of a list of case value pairs plus the variable index, e.g.,
        # [((x1, y1), 1), ((x3,y3,), 3), ...]
        # gaps are due to missing value screening
        # index can be used to retrieve the weight, if any, or variable properties.
        
        return self.dissimfunc(z, lenz)

    def euclid(self, z, lenz):
        """Compute and return Euclidean dissimilarity between two cases
        """

        if not self.weights:
            res = sum((item[0] - item[1])**2 for item, index in z)
        else:
            res = sum((item[0] - item[1])**2 * self.weights[index] for item, index in z)
        res = math.sqrt(res)
        if self.missing == Dissimfuncs.RESCALE:
            return res * (self.numclusvars / lenz)  # we know denom can't be zero
        else:  # just ignore missings
            return res

    def minkowski(self, z, lenz):
        """Compute and return Minkowski dissimilarity between two cases
        """

        if not self.weights:
            res = sum(abs(item[0] - item[1])**self.power for item, index in z)
        else:
            res = sum(abs(item[0] - item[1])**self.power * self.weights[index] for item, index in z)
        res = res**(1./self.power)

        if self.missing == Dissimfuncs.RESCALE:
            return res * (self.numclusvars / lenz)  # we know denom can't be zero
        else:  # just ignore missings
            return res
    
    def maxdiff(self, z, lenz):
        """maximum abs difference across all variables 
        If weighting, use weight of the first variable with maximum diff
        There is no rescaling."""
        
        themax = max(abs(item[0]-item[1]) for item, index in z)
        if not self.weights:
            themax = max(abs(item[0]-item[1]) for item, index in z)
            return themax
        else:
            themax = max((abs(item[0]-item[1]), index) for item, index in z)
            return themax[0] * self.weights[themax[1]]
    
    def absdiff(self, z, lenz):
        """sum of absolute differences"""
        
        if not self.weights:
            res = sum(abs(item[0]-item[1]) for item, index in z)
        else:
            res = sum(abs(item[0]-item[1])  * self.weights[index] for item, index in z)
        if self.missing == Dissimfuncs.RESCALE:
            return res * (self.numvars / lenz)
        else:
            return res

    def gower(self, z, lenz):
        """Gower index
        
        Assumes variables are prepared and vrange has been calculated with constant variables
        already failed.
        """
        
        cum = 0.
        cumwt = 0.
        
        for item, index in z:
            vl = self.varmeta[index][1]
            if vl != "nominal":
                dis = abs(item[0] - item[1])/self.vrange[index]
            else:
                dis = (item[0] != item[1])   # count if nominal
            if self.weights:
                cum += dis * self.weights[index]
                cumwt += self.weights[index]
            else:
                cum += dis
                cumwt += 1
                
        return cum / cumwt

    
class Silcomp(object):
    """Compute silhouette statistics"""
    
    def __init__(self, cases, variables, dissim, missing, renumberordinal, weights, minkowskipower, vardict):
        """cases is a sequence of sequences where the first element of each is the cluster number
        cluster numbers must be positive integers.  Gaps are allowed, but they will be treated as
        clusters with size 0.  (This is useful if the data are just a sample of the clustering results)
        Data need not be sorted by cluster.
        variables is a list of clustering variables
        dissim is the keyword value of the function to be used to calculate the distances"""
        
        self.cases = cases
        self.numcases = len(self.cases)
        self.numvars = len(self.cases[0])  # number of variables used for clustering plus cluster id
        self.numclus = int(max(item[0] for item in self.cases if item[0] is not None))  # clusters could be missing in the data
        if not (2 <= self.numclus < self.numcases):
            raise ValueError(_("""Silhouettes cannot be computed: either there is only one cluster
or every cluster is size 1"""))
        self.diffmeas = Dissimfuncs(cases, variables, dissim, missing, renumberordinal, weights, minkowskipower,
            vardict)
        self.disstats = self.numcases * [0]  # one row for each case and one col for each cluster
        
    def calcsil(self):
        """Calculate dissimilarity statistics
        
        results go in casesil: one row per case containing ownclus, neighborclus, and dissimilarity"""
        
        # for each case, compute average distance from own-cluster members and for
        # each other cluster.  For smallest foreign cluster, use average distance for
        #silhouette calculation
        
        for c, case in enumerate(self.cases):
            clusstatsdiffs = self.numclus * [0]  # accumulate sum of diff measures by cluster
            clusstatscounts = self.numclus * [0]  # corresponding counts
            for i in range(self.numcases):
                diff = self.diffmeas.calc(case, self.cases[i])
                if not diff is None:
                    try:     # cluster number could be missing
                        clusnum = int(self.cases[i][0]) -1  #cluster numbers start at 1
                        clusstatsdiffs[clusnum] += diff
                        clusstatscounts[clusnum]+= 1
                    except:
                        pass
            self.disstats[c] = [Silcomp.saferatio(diff, count) for diff, count in zip(clusstatsdiffs, clusstatscounts)]

        self.casesil =self.numcases * [(0,0,0)]
        for c in range(self.numcases):
            try:
                ownclus = int(self.cases[c][0]) - 1
                a = self.disstats[c][ownclus]
                self.disstats[c][ownclus] = 1e100   # exclude from min calc
                try:
                    b = min([diff for diff in self.disstats[c] if not diff is None])  # next best
                    nextbestdiff = (b - a) / max(a, b)
                    self.casesil[c] = (ownclus+1, self.disstats[c].index(b)+1, nextbestdiff)
                except:
                    self.casesil[c] = (ownclus+1, None, None)
            except:
                self.casesil[c] = (None, None, None)
                
        del(self.cases)
    
    @staticmethod
    def saferatio(diff, count):
        if count == 0:
            return None
        else:
            return diff/count
        

# main routine
def sil(cluster, variables, nextbest, silhouette, dissim="euclid",  missing="rescale", histogram=True, counts=False,
    orientation="horizontal", bar=True, renumberordinal=False, minkowskipower=None, weights=None):
    """Calculate cluster silhouette statistics and optionally plot them

    cluster is the variable holding the cluster number for each case
    variables lists the variables used for clustering
    dissim is the dissimilarity function to use
"""
    # user missing values are converted to sysmis
    
    # debugging
    # makes debug apply only to the current thread
    #try:
        #import wingdbstub
        #if wingdbstub.debugger != None:
            #import time
            #wingdbstub.debugger.StopDebug()
            #time.sleep(2)
            #wingdbstub.debugger.StartDebug()
        #import thread
        #wingdbstub.debugger.SetDebugThreads({thread.get_ident(): 1}, default_policy=0)
        ## for V19 use
        ##    ###SpssClient._heartBeat(False)
    #except:
        #pass
    
    #info = NonProcPivotTable("INFORMATION", tabletitle=_("Information"),
        #columnlabels=["Count"])

    if cluster.lower() in [v.lower() for v in variables]:
        raise ValueError(_("""The cluster number variable cannot be included as a clustering variable"""))
    
    allvars = [cluster] + variables
    vardict = spssaux.VariableDict(allvars)
    if len(vardict.variables) != len(allvars):
        raise ValueError(_("""There are duplicate or undefined variables specified.  Note that variable names are case sensitive"""))

    cases = list(spssdata.Spssdata(allvars, names=False).fetchall())
    for c in range(len(cases)):
        cases[c] = list(cases[c])
        
    silcomp = Silcomp(cases, variables, dissim, missing, renumberordinal, weights, minkowskipower, vardict)
    silcomp.calcsil()
    
    # create result variables
    inputdataset = spss.ActiveDataset()
    if inputdataset == "*":   #unnamed - assign a name
        inputdataset = "D" + str(random.uniform(0,1))
        spss.Submit("DATASET NAME %s" % inputdataset)

    with spss.DataStep():
        ds = spss.Dataset()
        nbindex = createvar(ds, nextbest, _("Next Best Cluster"), [5, 8, 0], "NOMINAL")
        silindex = createvar(ds, silhouette, _("Silhouette"), [5, 8, 3], "SCALE")
        
        for i, case in enumerate(silcomp.casesil):
            ds.cases[i, nbindex] = silcomp.casesil[i][1]   # next best cluster
            ds.cases[i, silindex] = silcomp.casesil[i][2]     # silhouette statistic
    
    # aggregate results and retrieve
    aggdsname = "D" + str(random.uniform(0,1))
    spss.Submit(r"""DATASET DECLARE %(aggdsname)s.
TEMPORARY.
SELECT IF not SYSMIS(%(silhouette)s).
AGGREGATE
  /OUTFILE='%(aggdsname)s'
  /BREAK=%(cluster)s
  /N_BREAK=N
  /silhouette_mean=MEAN(%(silhouette)s) 
  /silhouette_min=MIN(%(silhouette)s) 
  /silhouette_max=MAX(%(silhouette)s).""" % locals())
    
    aggstats = sorted(spssdata.Spssdata(dataset=aggdsname).fetchall())
    spss.Submit("DATASET CLOSE %s" % aggdsname)
    rowlabels = [str(int(item[0])) for item in aggstats if item[0] is not None]
    cells = [item[1:] for item in aggstats if item[0] is not None]
    try:
        tcases = sum(item[1] for item in aggstats if item[1] is not None)
        tsum = sum(item[2] * item[1] for item in aggstats if item[1] is not None)
        tmean = tsum / tcases
        tmin = min(item[3] for item in aggstats if item[3] is not None)
        tmax = max(item[4] for item in aggstats if item[4] is not None)
        rowlabels.append(_("""Total"""))
        cells.append([tcases, tmean, tmin, tmax])
    except:
        pass
    
    StartProcedure("Cluster Silhouettes", "STATSCLUSSIL")
    tbl = spss.BasePivotTable(_("""Silhouette Statistics"""), "STATSCLUSSILSUM", isSplit=False,
        caption= _("""Dissimilarity measure = %s""" % silcomp.diffmeas.dissimlabel))
    tbl.SimplePivotTable(rowdim = _("""Cluster"""), rowlabels = rowlabels,
        coldim = _("""Statistics"""),
        collabels = [_("""Case Count"""), _("""Mean"""), _("""Minimum"""), _("""Maximum""")],
        cells = cells)
    spss.EndProcedure()
    
    spss.Submit("DATASET ACTIVATE %s." % inputdataset)
    if (histogram or bar or counts):
        paneling = (orientation == "horizontal" and "across") or "down"
        label = _("""Histogram of Silhouettes""")
        footnote = _("""Dissimilarity measure = %s""") % silcomp.diffmeas.dissimlabel
        title = _("""Silhouettes by Cluster""")
        
        if histogram:
            spss.Submit(r"""GGRAPH
/GRAPHDATASET NAME="graphdataset"
VARIABLES=%(cluster)s [LEVEL=nominal] %(silhouette)s [LEVEL=ratio] 
MISSING=LISTWISE REPORTMISSING=NO
/GRAPHSPEC SOURCE=VIZTEMPLATE(NAME="Histogram"[LOCATION=LOCAL]
MAPPING( "x"="%(silhouette)s"[DATASET="graphdataset"] 
"Footnote"="%(footnote)s" "Title"="%(title)s"
"Panel "+
"%(paneling)s"="%(cluster)s" [DATASET="graphdataset"] "Summary"="count"))
VIZSTYLESHEET="Traditional"[LOCATION=LOCAL]
    LABEL='%(label)s'
    DEFAULTTEMPLATE=NO.""" % locals())
            
        if bar:
            title = _("""Mean Silhouette by Cluster""")
            spss.Submit(r"""GGRAPH
/GRAPHDATASET NAME="graphdataset"
  VARIABLES=%(cluster)s [LEVEL=nominal] %(nextbest)s [LEVEL=nominal] %(silhouette)s [LEVEL=ratio] 
  MISSING=LISTWISE REPORTMISSING=NO
/GRAPHSPEC SOURCE=VIZTEMPLATE(NAME="3-D Bar"[LOCATION=LOCAL]
  MAPPING( "x"="%(cluster)s" [DATASET="graphdataset"] "y"="%(silhouette)s" [DATASET="graphdataset"] 
  "z"="%(nextbest)s"[DATASET="graphdataset"] "Summary"="mean" "Footnote"='%(footnote)s' "Title"="%(title)s"))
  VIZSTYLESHEET="Traditional"[LOCATION=LOCAL]
  LABEL='3-D BAR: Silhouettes'
  DEFAULTTEMPLATE=NO.""" % locals())
            
        if counts:
            xlabel = _("""Assigned Cluster""")
            ylabel = _("""Percentage of Assigned Cluster""")
            zlabel = _("""Next Best Cluster""")
            title =_("""Distribution of Next Best Clusters""")
            spss.Submit(r"""GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=%(cluster)s COUNT()[name="COUNT"] %(nextbest)s 
    MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: %(cluster)s=col(source(s), name("%(cluster)s"), unit.category())
  DATA: COUNT=col(source(s), name("COUNT"))
  DATA: %(nextbest)s=col(source(s), name("%(nextbest)s"), unit.category())
  COORD: rect(dim(1,2,3))
  GUIDE: axis(dim(1), label("%(zlabel)s"))
  GUIDE: axis(dim(2), label("%(xlabel)s"))
  GUIDE: axis(dim(3), delta(20), label("%(ylabel)s"))
  GUIDE: text.title(label("%(title)s"))
  GUIDE: text.footnote(label("%(footnote)s"))
  SCALE: linear(dim(3), min(0), max(100))
  ELEMENT: interval(position(summary.percent(%(nextbest)s*%(cluster)s*COUNT, base.coordinate(dim(2)))), 
    shape.interior(shape.square))
END GPL.""" % locals())

    
def createvar(ds, name, varlabel, format, measlevel):
    """Create or modify a numeric variable in an open dataset and return index
    
    ds is the dataset object
    name, varlabel, format, and measlevel are the variable properties
    format must be a sequence of type, width, decimals"""
    
    try:
        index = ds.varlist[name].index
        ds.varlist[index].type = 0    # just in case existing variable was a string
    except:
        ds.varlist.append(name, 0)
        index = len(ds.varlist) - 1
    ds.varlist[index].label = varlabel
    ds.varlist[index].format = format
    ds.varlist[index].measurementLevel = measlevel.upper()
    return index


def Run(args):
    """Execute the STATS CLUSTER SIL command"""

    args = args[list(args.keys())[0]]
    ###print args   #debug
    
    ###debugging
    #try:
        #import wingdbstub
        #if wingdbstub.debugger != None:
            #import time
            #wingdbstub.debugger.StopDebug()
            #time.sleep(2)
            #wingdbstub.debugger.StartDebug()
    #except:
        #pass

    oobj = Syntax([
        Template("VARIABLES", subc="",  ktype="existingvarlist", var="variables", islist=True),
        Template("CLUSTER", subc="",  ktype="existingvarlist", var="cluster", islist=False),
        Template("DISSIMILARITY", subc="",  ktype="str", var="dissim", vallist=["euclid", "absdiff", "maximum",
            "gower", "minkowski"]),
        Template("MINKOWSKIPOWER", subc="", ktype="float", var="minkowskipower"),
        Template("NEXTBEST", subc="", ktype="varname", var="nextbest"),
        Template("SILHOUETTE", subc="", ktype="varname", var="silhouette"),
        Template("RENUMBERORDINAL", subc="OPTIONS", ktype="bool", var="renumberordinal"),
        Template("VARWEIGHTS", subc="OPTIONS", ktype="float", var="weights", islist=True, vallist=[0]),
        Template("MISSING", subc="OPTIONS", ktype="str", var="missing", vallist=['rescale','ignore', 'omit']),
        Template("HISTOGRAM", subc="OUTPUT", ktype="bool", var="histogram"),
        Template("ORIENTATION", subc="OUTPUT", ktype="str", var="orientation", vallist=["horizontal", "vertical"]),
        Template("THREEDBAR", subc="OUTPUT", ktype="bool", var="bar"),
        Template("THREEDCOUNTS", subc="OUTPUT", ktype="bool", var="counts"),
        Template("HELP", subc="", ktype="bool")])
    
        # ensure localization function is defined
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg

        # A HELP subcommand overrides all else
    if "HELP" in args:
        #print helptext
        helper()
    else:
            processcmd(oobj, args, sil, vardict=spssaux.VariableDict())

def helper():
    """open html help in default browser window
    
    The location is computed from the current module name"""
    
    import webbrowser, os.path
    
    path = os.path.splitext(__file__)[0]
    helpspec = "file://" + path + os.path.sep + \
         "markdown.html"
    
    # webbrowser.open seems not to work well
    browser = webbrowser.get()
    if not browser.open_new(helpspec):
        print(("Help file not found:" + helpspec))
try:    #override
    from extension import helper
except:
    pass        
#class NonProcPivotTable(object):
    #"""Accumulate an object that can be turned into a basic pivot table once a procedure state can be established"""
    
    #def __init__(self, omssubtype, outlinetitle="", tabletitle="", caption="", rowdim="", coldim="", columnlabels=[],
                 #procname="Messages"):
        #"""omssubtype is the OMS table subtype.
        #caption is the table caption.
        #tabletitle is the table title.
        #columnlabels is a sequence of column labels.
        #If columnlabels is empty, this is treated as a one-column table, and the rowlabels are used as the values with
        #the label column hidden
        
        #procname is the procedure name.  It must not be translated."""
        
        #attributesFromDict(locals())
        #self.rowlabels = []
        #self.columnvalues = []
        #self.rowcount = 0

    #def addrow(self, rowlabel=None, cvalues=None):
        #"""Append a row labelled rowlabel to the table and set value(s) from cvalues.
        
        #rowlabel is a label for the stub.
        #cvalues is a sequence of values with the same number of values are there are columns in the table."""
        
        #if cvalues is None:
            #cvalues = []
        #self.rowcount += 1
        #if rowlabel is None:
            #self.rowlabels.append(str(self.rowcount))
        #else:
            #self.rowlabels.append(rowlabel)
        #if not spssaux._isseq(cvalues):
            #cvalues = [cvalues]
        #self.columnvalues.extend(cvalues)
        
    #def generate(self):
        #"""Produce the table assuming that a procedure state is now in effect if it has any rows."""
        
        #privateproc = False
        #if self.rowcount > 0:
            #try:
                #table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            #except:
                #StartProcedure(_("Messages"), self.procname)
                #privateproc = True
                #table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            #if self.caption:
                #table.Caption(self.caption)
            ## Note: Unicode strings do not work as cell values in 18.0.1 and probably back to 16
            #if self.columnlabels != []:
                #table.SimplePivotTable(self.rowdim, self.rowlabels, self.coldim, self.columnlabels, self.columnvalues)
            #else:
                #table.Append(spss.Dimension.Place.row,"rowdim",hideName=True,hideLabels=True)
                #table.Append(spss.Dimension.Place.column,"coldim",hideName=True,hideLabels=True)
                #colcat = spss.CellText.String("Message")
                #for r in self.rowlabels:
                    #cellr = spss.CellText.String(r)
                    #table[(cellr, colcat)] = cellr
            #if privateproc:
                #spss.EndProcedure()
                
def attributesFromDict(d):
    """build self attributes from a dictionary d."""
    self = d.pop('self')
    for name, value in d.items():
        setattr(self, name, value)
        

    """Manage a log file"""
    
    def __init__(self, logfile):
        """logfile is the file name or None"""

        self.logfile = logfile
        if self. logfile:
            self.file = open(logfile, "w")
            self.starttime = time.time()
            self.file.write("%.2f %s Starting log\n" % (time.time() - self.starttime, time.asctime()))
            
    def __enter__(self):
        return self
    
    def write(self, text):
        if self.logfile:
            self.file.write("%.2f: %s\n" % (time.time() - self.starttime,  text))
            self.file.flush()
            
    def close(self):
        if self.logfile:
            self.write("Closing log")
            self.file.close()

def StartProcedure(procname, omsid):
    """Start a procedure
    
    procname is the name that will appear in the Viewer outline.  It may be translated
    omsid is the OMS procedure identifier and should not be translated.
    
    Statistics versions prior to 19 support only a single term used for both purposes.
    For those versions, the omsid will be use for the procedure name.
    
    While the spss.StartProcedure function accepts the one argument, this function
    requires both."""
    
    try:
        spss.StartProcedure(procname, omsid)
    except TypeError:  #older version
        spss.StartProcedure(omsid)