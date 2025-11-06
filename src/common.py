import streamlit as st
import pandas as pd
import io
import uuid
import os

def clear_cache_button():
   if st.button("Clear Cache"):
        # Clear cache for both newer and older Streamlit versions
        if hasattr(st, "cache_data"):
            st.cache_data.clear()
        if hasattr(st, "cache_resource"):
            st.cache_resource.clear()
        st.success("Cache cleared!")

    # initialize global session state variables if not already present
    # DataFrames

def v_space(n, col=None):
    for _ in range(n):
        if col:
            col.write("")
        else:
            st.write("")

def page_setup():
    # streamlit configs
    st.set_page_config(
        page_title="Corromics",
        page_icon="assets/corromics_icon.png",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )
    for key in dataframe_names:
        if key not in st.session_state:
            st.session_state[key] = pd.DataFrame()
    if "data_preparation_done" not in st.session_state:
        st.session_state["data_preparation_done"] = False

    with st.sidebar:
        with st.expander("⚙️ Settings", expanded=True):
            st.selectbox(
                "image export format",
                ["svg", "png", "jpeg", "webp"],
                key="image_format",
            )
        v_space(1)
        # Add the clear cache button
        clear_cache_button()
        v_space(1)
        
        try:
            st.image("https://raw.githubusercontent.com/abzer005/Corromics/main/assets/corromics_full_logo.png",
                     use_container_width=True)
        except TypeError:
            st.image("https://raw.githubusercontent.com/abzer005/Corromics/main/assets/corromics_full_logo.png", 
                     use_column_width=True)

dataframe_names = ("md",
                   "ft",
                   "nw",
                   "an_gnps",
                   "an_analog")


def reset_dataframes():
    for key in dataframe_names:
        st.session_state[key] = pd.DataFrame()


def open_df(file):
    separators = {"txt": "\t", "tsv": "\t", "csv": ","}
    try:
        if type(file) == str:
            ext = file.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        else:
            ext = file.name.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        
        # sometimes dataframes get saved with unnamed index, that needs to be removed
        if "Unnamed: 0" in df.columns:
            df.drop("Unnamed: 0", inplace=True, axis=1)
        return df
    except:
        return pd.DataFrame()
    
import pandas as pd

def open_stats_ft(file):
    """
    Read FBMN-Stats feature table without dropping the first (index) column.
    Unlike open_df(), this keeps 'Unnamed: 0' so we can use it as sample IDs.
    """
    separators = {"txt": "\t", "tsv": "\t", "csv": ","}
    try:
        if isinstance(file, str):
            ext = file.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators.get(ext, ","))
            else:
                df = pd.read_excel(file)
        else:
            ext = file.name.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators.get(ext, ","))
            else:
                df = pd.read_excel(file)

        return df
    except Exception:
        return pd.DataFrame()


def show_table(df, title="", col="", download=True):
    if col:
        col = col
    else:
        col = st
    if download:
        col.download_button(
            f"Download Table",
            df.to_csv(sep="\t").encode("utf-8"),
            title.replace(" ", "-") + ".tsv",
            key=uuid.uuid1(),
        )
    col.dataframe(df, use_container_width=True)


def show_fig(fig, download_name, container_width=True):

    # Set default image format to 'svg' if not specified in session state
    image_format = st.session_state.get('image_format', 'svg')

    st.plotly_chart(
        fig,
        use_container_width=container_width,
        config={
            "displaylogo": False,
            "modeBarButtonsToRemove": [
                "zoom",
                "pan",
                "select",
                "lasso",
                "zoomin",
                "autoscale",
                "zoomout",
                "resetscale",
            ],
            "toImageButtonOptions": {
                "filename": download_name,
                "format": image_format,
            },
        },
    )


def download_plotly_figure(fig, filename="", col=""):
    buffer = io.BytesIO()
    fig.write_image(file=buffer, format="png")

    if col:
        col.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )
    else:
        st.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )

  
def get_max_correlation_limit():
    limits = {
        "restricted": 1_000_000, # Safe limit for restricted (cloud/GNPS/streamlit)
        "local": 10_000_000,  #Safe upper limit for local mode
    }
    return limits["restricted"] if st.session_state.is_restricted_mode else limits["local"]


def initialize_app():

    # Try to read from Streamlit secrets
    environment_mode = st.secrets.get("ENVIRONMENT_MODE", None)

    # If not found, fallback to environment variable
    if environment_mode is None:
        environment_mode = os.getenv("ENVIRONMENT_MODE", "local")
    
    # Store in session state for global access
    st.session_state.is_restricted_mode = environment_mode in ["gnps", "streamlit"]

    # Optional: warning if restricted mode is active
    with st.sidebar:
        if st.session_state.is_restricted_mode:
            st.warning(f"⚠️ Restricted mode is active: \nMax correlation limit: **{get_max_correlation_limit():,}**")
        else:
            st.success(f"App running in **local mode**.\nMax correlation limit: **{get_max_correlation_limit():,}**")
