#!/usr/bin/env python


from __future__ import print_function
import os
import sys
import argparse
parser = argparse.ArgumentParser('''
This script runs the following checks on the specified datacard:
- impact computation
- significance estimation for discovery (with toys)
- significance estimation for discovery (with asymptitic formulae)
- upper limits
- log-likelihood profiles vs. parameters
Details of the commands in the code.''')
parser.add_argument('-d', '--datacard', required=True           , help='path to the combine datacard to open' )
parser.add_argument('-p', '--pranges' , default=''              , help='parameter ranges to use in the form "par1=min,max:par2=min,max"')
parser.add_argument('-e', '--expected', default='1'             , help='expected signal strength for the Asimov generation')
parser.add_argument('-l', '--label'   , default='combine_checks', help='used to label the output directory')
parser.add_argument('-r', '--rmin'    , default='0'             , help='lower boundary for r')
parser.add_argument('-R', '--rmax'    , default='10'            , help='higher boundary for r')
parser.add_argument('-t', '--toys'    , default='10000'         , help='toys to use to compute the significance or the upper limit')
parser.add_argument('-U', '--unblind' , action='store_true'     , help='run unblinded tests')
parser.add_argument(      '--log'     , action='store_true'     , help='store stdout and stderr in log files instead of printing them')
parser.add_argument('-I', '--impacts' , action='store_true'     , help='run the impacts computation')
parser.add_argument('-A', '--sig-asym', action='store_true'     , help='run the significance computation using asymptotic formulae')
parser.add_argument('-T', '--sig-toys', action='store_true'     , help='run the significance computation using toys')
parser.add_argument('-L', '--lhscan'  , action='store_true'     , help='run the likelihood scan')
parser.add_argument('-D', '--teststat', action='store_true'     , help='compute the test statistics distribution')
parser.add_argument('-F', '--fit'     , action='store_true'     , help='compute the best fit for the signal strength')
parser.add_argument('-G', '--gridscan', action='store_true'     , help='scan grid points for later limit computation')
parser.add_argument('-P', '--rpoint'  , default=None            , help='r point to be used with --gridscan')
parser.add_argument(      '--lhparams', default='r'             , help='parameters to run the likelihood scan in the form "par1 par2 par3"')
parser.add_argument(      '--limit'   , action='store_true'     , help='run the upper limit with toys')
parser.add_argument(      '--alimit'  , action='store_true'     , help='run the upper limit with asymptotic formulae')
parser.add_argument(      '--cl'      , default='0.9'           , help='confidence level used for the upper limit')
parser.add_argument(      '--grid'    , default='0.5'           , help='quantile considered for the upper limit computation')
parser.add_argument(      '--tofreeze', default=None            , help='parameters to freeze to help the fit convergence')
parser.add_argument('-m', '--method'  , default='CLs'           , help='method to use for upper limits', choices=['CLs', 'CLsplusb'])
parser.add_argument(      '--rebin'   , default=None            , help='rebin the mass distribution when running FitDiagnostic fits')
parser.add_argument('--generate-nuisances', action='store_true' , help='equivalent of combine argument --generateNuisances')
args = parser.parse_args()

class Command:
  WARNING = lambda self, msg: print("\033[93m[WARNING] "+ msg + "\033[0m")
  ERROR   = lambda self, msg: print("\033[91m[ERROR] "  + msg + "\033[0m")
  INFO    = lambda self, msg: print("\033[92m[INFO]  "  + msg + "\033[0m")
  CODE    = lambda self, msg: print("\033[95m" + msg + "\033[0m")
  def __init__(self, name, commands, output):
    assert type(commands)==list, "Please pass command(s) via a python list"
    self.n = name
    self.o = output
    self.c = commands
    self.w = os.getcwd()
    self.f = False  #failed
  
  def run(self, runme=True):
    if not runme: return

    self.INFO("running task {}".format(self.n))
    self.INFO("the output directory is")
    print(self.o)

    if not os.path.exists(self.o):
      os.makedirs(self.o)

    os.chdir(self.o)
    for i, cmd in enumerate(self.c):
      cmd = cmd+'>> {O}/LOGFILE.txt'.format(O=self.o) if args.log else cmd
      print("   cmd -====  ", cmd)
      self.CODE(">> "+cmd)
      if args.log:
        self.INFO("logging stdout and stderr into")
        print(self.o+'/LOGFILE.txt')

      if args.log:
        os.system('echo \>\> "%s 2>&1"' %cmd)
      ret = os.system(cmd)
      if ret:
        self.ERROR("an error was encountered when running the last command")
        self.f = True
    
    if self.f:
      self.WARNING("the {} task encountered some errors".format(self.n))
    else:
      self.INFO("task {} run succesfully".format(self.n))
    os.chdir(self.w)

DATACARD   = os.path.abspath(args.datacard)
WORKSPACE  = os.path.basename(args.datacard).replace('.txt', '.root')
PARAMETERS = '--setParameterRanges "{}"'.format(args.pranges) if args.pranges!='' else ''
PARAMETERS = PARAMETERS + ' --freezeParameters "{}"'.format(args.tofreeze) if args.tofreeze is not None else PARAMETERS
BLINDER    = '-t -1 --expectSignal {EXP}'.format(EXP=args.expected) if not args.unblind else ''
MODE       = "--generateNuisances=1 --generateExternalMeasurements=0 --fitNuisances=0 --testStat LHC" if args.generate_nuisances else '--LHCmode LHC-limits'
print(' TEst')
OUTPUT = os.path.abspath('_'.join([args.label, 'rMin'+args.rmin, 'rMax'+args.rmax]+['unblinded' if args.unblind else 'rAsimov'+args.expected]))
OUTPUT_IMPACTS  = OUTPUT+'/impacts'
OUTPUT_SIG_TOYS = OUTPUT+'/significance_toys'
OUTPUT_SIG_ASYM = OUTPUT+'/significance_asymptotic'
OUTPUT_LHSCAN   = OUTPUT+'/lh_scan'
OUTPUT_TESTSTAT = OUTPUT+'/teststat'
OUTPUT_LIMIT    = '_'.join([OUTPUT+'/limit', 'CL'+args.cl, 'grid'+args.grid.replace('.','p')])
OUTPUT_ALIMIT   = '_'.join([OUTPUT+'/asymptoticlimit', 'CL'+args.cl])
OUTPUT_FIT      = OUTPUT+'/bestfit'
OUTPUT_GRIDSCAN = OUTPUT+'/gridscan{}'.format(args.rpoint.replace('.','p')) if args.rpoint is not None else None



CMD_IMPACTS = '\n'.join([
  'combine -M FitDiagnostics -d {DAT} {B} --plots --rMin "{r}" --rMax "{R}" --minos all {PAR} {REB}',
  'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnosticsTest.root -g plots.root',
  'text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName data_obs --X-assign-flatParam-prior',
  'combineTool.py -M Impacts -d {WSP} -m 1.777 {B} --doInitialFit --robustFit=1 --rMin "{r}" --rMax "{R}"',
  'combineTool.py -M Impacts -d {WSP} -m 1.777 {B} --doFits --parallel 5 --rMin "{r}" --rMax "{R}"',
  'combineTool.py -M Impacts -d {WSP} -m 1.777 {B} --rMin "{r}" --rMax "{R}" -o impacts.json',
  'plotImpacts.py -i impacts.json -o impacts',
]).format(
  DAT=DATACARD      ,
  EXP=args.expected ,
  r  =args.rmin     ,
  R  =args.rmax     ,
  RAN=args.pranges  ,
  WSP=WORKSPACE     ,
  PAR=PARAMETERS    ,
  B  =BLINDER       ,
  REB='--rebinFactor {}'.format(args.rebin) if args.rebin is not None else '',
).split('\n')

#print("Impacts command  ", CMD_IMPACTS)

CMD_SIGNIFICANCE_TOYS = '\n'.join([
  'combine -M HybridNew -d {DAT} --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult  -T {T} {PAR}',
]).format(
  T  =args.toys ,
  DAT=DATACARD  ,
  PAR=PARAMETERS,
).split('\n')

CMD_SIGNIFICANCE_ASYM = '\n'.join([
  'combine -M Significance -d {DAT}',
]).format(
  DAT=DATACARD
).split('\n')

# then plot with
# plot1DScan.py /path/to/main.root --others /path/to/additional.root --POI r --output myscan --main-label main_file
CMD_LHSCAN = '\n'.join([
  'text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName data_obs',
  'for nui in {NUI}; do \
    echo scanning $nui;\
    combine {WSP} -M MultiDimFit --algo grid --points 20 {PAR} -m 1.777 --rMin "{r}" --rMax "{R}" -P $nui -n scan_$nui {BLI} --squareDistPoiStep --autoRange 3;\
    plot1DScan.py higgsCombinescan_$nui.MultiDimFit.mH1.777.root --POI  \"$nui\" --output plot_$nui;\
  done'
]).format(
  DAT=DATACARD  ,
  WSP=WORKSPACE ,
  BLI=BLINDER   ,
  PAR=PARAMETERS,
  r  =args.rmin ,
  R  =args.rmax ,
  NUI=args.lhparams,
).split('\n')

CMD_TESTSTAT = '\n'.join([
  'text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName data_obs',
  'for n in {{1..100}}; do \
    combine -M HybridNew {WSP} --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult -T {TOY} {PAR} --rMin {r} --rMax {R} -s -1; \
  done',
  'hadd -k -f grid.root *.root',
  '$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotTestStatCLs.py grid.root {BLI} --val all --mass 120',
]).format(
  DAT=DATACARD ,
  WSP=WORKSPACE,
  TOY=args.toys,
  r  =args.rmin,
  R  =args.rmax,
  PAR=PARAMETERS,
  BLI='-E -q {}'.format(args.grid) if not args.unblind else '',
).split('\n')

CMD_LIMIT = '\n'.join([
  'text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName data_obs --X-assign-flatParam-prior',
#  'text2workspace.py {DAT} -m 1.777 -o {WSP} -D data_obs',
#  "combine -M HybridNew {WSP} {BLI} {MOD} --X-rtd MINIMIZER_freezeDisassociatedParams -T {TOY} -C {CL}  --plot='limit_combined_hybridnew_{CL}.pdf' --rMin {r} --rMax {R} {PAR} --rule {RL}"
#  "combine -M HybridNew {WSP} {BLI} {MOD} --X-rtd MINIMIZER_freezeDisassociatedParams -T {TOY} -C {CL}  --plot='limit_combined_hybridnew_{CL}.pdf' --rMin {r} --rMax {R} {PAR} --rule {RL}"

  "combine -M HybridNew {WSP} {BLI} {MOD} --saveHybridResult  -T {TOY} -C {CL}  --plot='limit_combined_hybridnew_{CL}.pdf' --rMin {r} --rMax {R} {PAR} --rule {RL}"
#  "combine -M HybridNew {WSP} {BLI}  --generateNuisances={GEN} --generateExternalMeasurements=0 --fitNuisances=0 --testStat LHC --X-rtd MINIMIZER_freezeDisassociatedParams -T {TOY} -C {CL}  --plot='limit_combined_hybridnew_{CL}.pdf' --rMin {r} --rMax {R} {PAR} --rule {RL}"
]).format(
 #  GEN='1' if args.generate_nuisances else '0',
  DAT=DATACARD ,
  WSP=WORKSPACE,
  TOY=args.toys,
  CI =args.grid,
  CL =args.cl  ,
  RL =args.method,
  r  =args.rmin,
  R  =args.rmax,
  PAR=PARAMETERS,
  BLI='--expectedFromGrid {CI}'.format(CI=args.grid) if not args.unblind else '',
  MOD=MODE,
).split('\n')

CMD_FIT = '\n'.join([
  'text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName data_obs',
  'combine -M FitDiagnostics {WSP} --plots --rMin {r} --rMax {R} {PAR}',
]).format(
  WSP=WORKSPACE,
  DAT=DATACARD ,
  r  =args.rmin,
  R  =args.rmax,
  PAR=PARAMETERS,
).split('\n')

CMD_ALIMIT = '\n'.join([
  'combine -M AsymptoticLimits --run {BLI} --cl {CL} -d {DAT}'
]).format(
  BLI='blind' if not args.unblind else 'observed',
  CL =args.cl ,
  DAT=DATACARD,
).split('\n')

CMD_GRIDSCAN= '\n'.join([
  "text2workspace.py {DAT} -m 1.777 -o {WSP} --dataMapName=data_obs {FLT}",
  "combine {WSP} -M HybridNew -m 1.777 {MOD} --singlePoint {X} --saveToys --saveHybridResult -T {TOY} --clsAcc 0 --rMin {r} --expectedFromGrid 0.5 --rMax {R} {PAR}",
#  "combine {WSP} -M HybridNew -m 1.777 --LHCmode LHC-limits --singlePoint {X} --saveToys --saveHybridResult -T {TOY} --clsAcc 0 --rMin {r} --expectedFromGrid 0.5 --rMax {R} {PAR} --rule {RL}",
]).format(
  DAT=DATACARD  ,
  WSP=WORKSPACE ,
  GEN='1' if args.generate_nuisances else '0',
  TOY=args.toys,
  R  =args.rmax,
  r  =args.rmin,
  PAR=PARAMETERS,
  RL =args.method,
  X  =args.rpoint,
  MOD=MODE,
  FLT="--X-assign-flatParam-prior" if args.generate_nuisances else ""
).split('\n')

impacts_cmd  = Command('Impacts'                , CMD_IMPACTS           , OUTPUT_IMPACTS ).run(args.impacts )
print("-=----- ", impacts_cmd)
sig_asym_cmd = Command('Asymptotic significance', CMD_SIGNIFICANCE_ASYM , OUTPUT_SIG_ASYM).run(args.sig_asym)
sig_toys_cmd = Command('Toys significance'      , CMD_SIGNIFICANCE_TOYS , OUTPUT_SIG_TOYS).run(args.sig_toys)
lhscan_cmd   = Command('LHScan'                 , CMD_LHSCAN            , OUTPUT_LHSCAN  ).run(args.lhscan  )
limit_cmd    = Command('Upper limit'            , CMD_LIMIT             , OUTPUT_LIMIT   ).run(args.limit   )
alimit_cmd   = Command('Upper limit'            , CMD_ALIMIT            , OUTPUT_ALIMIT  ).run(args.alimit  )
teststat_cmd = Command('Test-stat distribution' , CMD_TESTSTAT          , OUTPUT_TESTSTAT).run(args.teststat)
fit_cmd      = Command('Test-stat distribution' , CMD_FIT               , OUTPUT_FIT     ).run(args.fit     )
grid_cmd     = Command('grid scan'              , CMD_GRIDSCAN          , OUTPUT_GRIDSCAN).run(args.gridscan)
