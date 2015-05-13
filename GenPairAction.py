import BareSquarer
import DavidSquarer
import IlkkaSquarer

def run(pa_objects):
    for pa_object in pa_objects:
        if pa_object['squarer']['type'] is 'Ilkka':
            print '\n**** Performing squaring ****\n'
            print 'WARNING: Assuming 3D Coulomb interaction'
            IlkkaSquarer.Square(pa_object)
            print '\n**** Performing breakup ****\n'
            IlkkaSquarer.Breakup(pa_object)
        elif pa_object['squarer']['type'] is 'David':
            print '\n**** Performing breakup ****\n'
            DavidSquarer.Breakup(pa_object)
            print '\n**** Performing squaring ****\n'
            DavidSquarer.Square(pa_object)
        elif pa_object['squarer']['type'] is 'None':
            print '\n**** Performing breakup ****\n'
            BareSquarer.Breakup(pa_object)
