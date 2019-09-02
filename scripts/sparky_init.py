def initialize_session(session):
    def vv_command(s=session):
        import paramagpy_sparky_macro
        paramagpy_sparky_macro.read_write_pcs_files(s)

    session.add_command("vv", "Read and Write PCS Files", vv_command)
