# net growth of seagrass accounting for temperature and light effects on photosynthesis and temperature dependent respiration
Net.growth.shade <- function(T., I., PT.max, PL.max, T.opt, T.max, Ik, R.max, RT.opt, RT.max, Ca, B.max){
  photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca, B.max) - Resp.temp(R.max, T., RT.opt, RT.max)
}
