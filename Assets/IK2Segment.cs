using UnityEngine;

namespace InverseKinematics2D
{
    [ExecuteInEditMode]
    public class IK2Segment : MonoBehaviour
    {
        public float length;
        public Interval angleLimit;
        [Range(0, 1)] public float targetAngle = 0.5f;

        public bool debug = true;

        [Header("Lizard Extra")]
        [Range(0.0f, 1.0f)] public float tension = 0.1f;

        private Color debugColor;

        public float debugTension;

        private void OnValidate()
        {
            debugColor = new Color(Random.Range(0.0f, 1.0f), Random.Range(0.0f, 1.0f), Random.Range(0.0f, 1.0f), 1.0f);
        }

        public void Init(float length, float angleMin, float angleMax)
        {
            this.length = length;
            this.angleLimit = new Interval(angleMin, angleMax);
        }

        private void Update()
        {
            if (transform.childCount > 0)
            {
                transform.GetChild(0).transform.localPosition = new Vector3(length, 0, 0);
            }
            debugTension = AppliedTension();
        }

        public float AppliedTension()
        {
            float angle = Vector2.SignedAngle(transform.parent.right, transform.right);
            float minT = MinTEdgeAngle;
            float maxT = MaxTEdgeAngle;

            float dmin = minT - angle;
            float dmax = angle - maxT;

            // outside tension zone
            if (dmin < 0 && dmax < 0)
            {
                return 0;
            }
            // on max side
            else if (dmin < 0)
            {
                return Mathf.Pow(dmax / (angleLimit.max - maxT), 4);
            }
            // on min side
            else
            {
                return Mathf.Pow(dmin / (minT - angleLimit.min), 4);
            }
        }

        private void OnDrawGizmos()
        {
            if (!debug)
            {
                return;
            }
            Color originalColor = Gizmos.color;
            Gizmos.color = debugColor;
            Gizmos.matrix = Matrix4x4.TRS(transform.position + transform.right * length / 2, transform.rotation, new Vector3(length, 0.2f, 0.2f));
            Gizmos.DrawCube(Vector3.zero, Vector3.one);
            Gizmos.matrix = Matrix4x4.identity;
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(angleLimit.min, transform.forward) * transform.parent.right * length);
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(angleLimit.max, transform.forward) * transform.parent.right * length);
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(MinTEdgeAngle, transform.forward) * transform.parent.right * length);
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(MaxTEdgeAngle, transform.forward) * transform.parent.right * length);
            Gizmos.color = originalColor;


            IK2 ik = GetComponentInParent<IK2>();
            ik?.DebugSegmentReach(this, ReachDebugMode.Back, debugColor);
        }

        public float HalfAngle => (angleLimit.min + angleLimit.max) / 2 - angleLimit.min;
        public float IdealAngle => (angleLimit.min + angleLimit.max) / 2;
        public float MinTEdgeAngle => IdealAngle - HalfAngle * (1 - tension);
        public float MaxTEdgeAngle => IdealAngle + HalfAngle * (1 - tension);
    }
}
