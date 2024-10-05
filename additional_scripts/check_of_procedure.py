# Verify the existence of the requested procedure
def check_of_procedure(procedure_name):
    return procedure_name in [
        "transcribe",
        "reverse",
        "complement",
        "reverse_complement",
    ]
