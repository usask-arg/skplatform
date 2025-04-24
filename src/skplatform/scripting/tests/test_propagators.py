import logging 

def test_propagator_variables():
    
    from skplatform.scripting.interface import default_tle
    from skplatform.scripting.instruments import SimulatorInstrument, AOSSky, TICFIRE
    
    for inst_class in [SimulatorInstrument, AOSSky, TICFIRE]:
        inst = inst_class(default_tle())
        inst.propagate()
        
        variables = inst.propagator._dataset_variables()
        ds = inst.get_dataset()
        # we don't care about the difference between coords and variables for this
        keys = [k for k in ds.keys()]
        keys.extend([c for c in ds.coords])
        for key in keys:
            if key in variables:
                variables = [x for x in variables if x != key]
            else:
                logging.error(f'Variable {inst_class.__name__}.{key} is defined in output but not in _dataset_variables()')
                raise KeyError
        
        if variables != []:
            logging.error(f'{inst_class.__name__}_dataset_variables() contains variables not exported in get_dataset(): {variables}')
            raise KeyError

if __name__ == "__main__":
    test_propagator_variables()