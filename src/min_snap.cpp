#include <min_snap/min_snap.h>


int main(int argc, const char** argv) {
    ros::init(argc, argv, "min_snap");
    ros::NodeHandle nh;
    min_snap ms(nh);
    ros::spin();
    return 0;
}