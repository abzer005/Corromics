import sys
import streamlit.web.cli as stcli

if __name__ == "__main__":
    sys.argv = ["streamlit", "run", "Home.py"]
    sys.exit(stcli.main())


# from streamlit.web import cli

# if __name__ == "__main__":
#     cli._main_run_clExplicit(
#         file="Home.py", command_line="streamlit run", args=[]
#     )
#     # we will create this function inside our streamlit framework