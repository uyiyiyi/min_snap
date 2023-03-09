#include <min_snap/min_snap.h>


int main(int argc, char** argv) {
    ros::init(argc, argv, "min_snap");
    ros::NodeHandle nh;
    min_snap ms(nh);    
    ms.visTraj(ms.Solution_x, ms.Solution_y, ms.Solution_z);
    // ros::spin();
    return 0;
}