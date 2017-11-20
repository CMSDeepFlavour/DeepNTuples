import pdb
import ROOT

from DataFormats.FWLite import Events, Handle

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.gSystem.Load("libDataFormatsPatCandidates.so")

print("Get Handles")

handle_trigger = Handle('edm::TriggerResults')
label_trigger = ("TriggerResults")

events = Events("Data/output1_0_1.root")

print("loop over events")

for iEvent in events:
    iEvent.getByLabel(label_trigger, handle_trigger)
    triggers = handle_trigger.product()
    triggerNames = triggers.triggerNames(handle_trigger)

    pdb.set_trace()

    for itrig in range(0, handle_trigger.size()):
        trigname = triggerNames.triggerName(itrig)
        print("trigname", trigname)

