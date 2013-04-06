import numpy
import binding


contacts_ts = numpy.array([0,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,1])
n = binding._count_binding_events_state_machine(contacts_ts)
print "Total number of binding events:", n

contacts_ts = numpy.array([0,0,0,1,1,1,0,0,0,1,1,1,0,2,1,4,5,3])
n = binding._count_binding_events_state_machine(contacts_ts)
print "Total number of binding events:", n
