from pysatAMISR.instruments import methods  # noqa F401

__all__ = ['isr_pf']

for inst in __all__:
    exec("from pysatAMISR.instruments import {:}".format(inst))
