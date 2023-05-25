using UnityEngine;

namespace InverseKinematics2D
{
    [ExecuteInEditMode]
    public class IK2Segment : MonoBehaviour
    {
        public float length;
        public Interval angleLimit;
        [Range(0, 1)] public float targetAngle = 0.5f;

        private Color debugColor;

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
        }

        private void OnDrawGizmos()
        {
            Color originalColor = Gizmos.color;
            Gizmos.color = debugColor;
            Gizmos.matrix = Matrix4x4.TRS(transform.position + transform.right * length / 2, transform.rotation, new Vector3(length, 0.2f, 0.2f));
            Gizmos.DrawCube(Vector3.zero, Vector3.one);
            Gizmos.matrix = Matrix4x4.identity;
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(angleLimit.min, transform.forward) * transform.parent.right);
            Gizmos.DrawRay(transform.position, Quaternion.AngleAxis(angleLimit.max, transform.forward) * transform.parent.right);
            Gizmos.color = originalColor;
        }
    }
}
