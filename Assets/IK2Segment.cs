using UnityEngine;

namespace InverseKinematics2D
{
    public class IK2Segment : MonoBehaviour
    {
        public float length;
        public Interval angleLimit;
        [Range(0, 1)] public float targetAngle = 0.5f;

        public void Init(float length, float angleMin, float angleMax)
        {
            this.length = length;
            this.angleLimit = new Interval(angleMin, angleMax);
        }
    }
}
