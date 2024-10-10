#
# Source this file to set up the environment to work with FIPPS_Softs
#

if [ -z "${PATH}" ]; then
   PATH=@FIPPS_Softs_BIN_DIR@; export PATH
else
   PATH=@FIPPS_Softs_BIN_DIR@:$PATH; export PATH
fi

if [ -n "${LD_LIBRARY_PATH}" ];  then
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:@FIPPS_Softs_LIB_DIR@  ; export LD_LIBRARY_PATH  # Linux, ELF HP-UX
fi

if [ -n "${DYLD_LIBRARY_PATH}" ];  then
  DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:@FIPPS_Softs_LIB_DIR@  ; export DYLD_LIBRARY_PATH  # Mac OS X
fi

if [ -n "${SHLIB_PATH}" ];  then
  SHLIB_PATH=$SHLIB_PATH:@FIPPS_Softs_LIB_DIR@  ; export SHLIB_PATH  # legacy HP-UX
fi

if [ -n "${LIBPATH}" ];  then
  LIBPATH=$LIBPATH:@FIPPS_Softs_LIB_DIR@  ; export LIBPATH  # AIX
fi

FIPPS_Softs_SYS=@FIPPS_Softs_INSTALL_DIR@ ; export FIPPS_Softs_SYS

source @FIPPS_Softs_BIN_DIR@/FIPPS_Softs-autocompletion.sh

echo ' ---> You are working now with FIPPS_Softs installed in ' @FIPPS_Softs_INSTALL_DIR@
