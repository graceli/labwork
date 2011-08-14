##
## Example script showing a way to add user-drawn geometry that updates
## as a trajectory is animated.
##
## To use this script:
## 1) load your trajectory
## 2) source frameupdate.vmd
## 3) enabletrace
## 4) do your thing :-)
## 5) disabletrace

proc enabletrace {} {
  global vmd_frame;
  trace variable vmd_frame([molinfo top]) w drawcounter
}

proc disabletrace {} {
  global vmd_frame;
  trace vdelete vmd_frame([molinfo top]) w drawcounter
}

proc drawcounter { name element op } {
  global vmd_frame;

  draw delete all
  # puts "callback!"
  draw color white
  set psperframe 1.03
  set time [format "%8.3f ps" [expr $vmd_frame([molinfo top]) * $psperframe]]
  draw text {10 10 0} "$time"
} 



