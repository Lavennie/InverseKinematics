using UnityEditor;
using UnityEngine;

namespace InverseKinematics1
{
    [ExecuteInEditMode]
    public class IK : MonoBehaviour
    {
        private readonly float ERROR = Mathf.Epsilon;

        public IKBone[] bones;
        public Transform endEffector;
        [Range(0.0f, 1.0f)] public float weight = 0.5f;

        private float fullLength;
        private float fullLengthSquared;

        //public Vector3 circlePoint;

        public void OnEnable()
        {
            fullLength = 0.0f;
            for (int i = 0; i < bones.Length; i++)
            {
                fullLength += bones[i].length;
            }
            fullLengthSquared = fullLength * fullLength;
            //circlePoint = Quaternion.Euler(Random.Range(-180.0f, 180.0f), Random.Range(-180.0f, 180.0f), Random.Range(-180.0f, 180.0f)) * (Vector2.right * Random.Range(0.0f, 1.0f));
        }

        public void Update()
        {
            //Debug.DrawLine(bones[0].transform.position, bones[0].transform.position + (bones[0].transform.rotation * (Vector2)circlePoint + bones[0].transform.forward) * 0.1f, Color.magenta, 0.1f);
            Resolve2Chain(0);

            //Vector3 a = new Vector3(bones[0].xLimit.min, 0, bones[0].zLimit.min);
            //Vector3 b = new Vector3(0, bones[0].yLimit.min, bones[0].zLimit.max);
            //bones[0].transform.rotation = Quaternion.Euler(Vector3.Lerp(a, b, weight));
            //bones[0].transform.rotation = bones[0].GetResting();
            //bones[1].transform.rotation = bones[0].transform.rotation * bones[1].GetResting();
            //bones[0].transform.rotation = Quaternion.Slerp(Quaternion.Euler(30, 0, 0), Quaternion.Euler(0, -40, 0), weight);
        }

        private void Resolve2Chain(int i)
        {
            // clamp end position to local yz plane
            Vector3 end = transform.parent.InverseTransformPoint(endEffector.position);
            end.x = 0;
            end = transform.parent.TransformPoint(end);
            float startAngle = Vector3.SignedAngle(transform.parent.forward, end - transform.parent.position, transform.parent.right);
            Debug.DrawLine(transform.position, end, Color.white, 0.01f);

            Vector3 dirToEnd = end - bones[i].transform.position;
            float distToEndSquared = dirToEnd.sqrMagnitude;
            float distToEnd = Mathf.Sqrt(distToEndSquared);
            dirToEnd /= distToEnd;
            if (distToEnd > bones[i].length + bones[i + 1].length)
            {
                bones[i].transform.forward = dirToEnd;
                bones[i + 1].transform.forward = dirToEnd;
            }
            else if (distToEnd < Mathf.Abs(bones[i].length - bones[i + 1].length))
            {
                bones[i].transform.forward = dirToEnd;
                bones[i + 1].transform.forward = -dirToEnd;
            }
            else
            {
                // law of cosines
                float angle1 = Mathf.Rad2Deg * Mathf.Acos((bones[i].length * bones[i].length + distToEnd * distToEnd -
                    bones[i + 1].length * bones[i + 1].length) / (2.0f * bones[i].length * distToEnd));
                AngleRegion region1 = new AngleRegion(bones[i].xLimit);
                AngleRegion angleRegion1 = new AngleRegion(-angle1, angle1);
                region1.Intersect(angleRegion1);
                angle1 = Mathf.Lerp(region1.min, region1.max, weight);
                DebugRegion(i - 1, region1, Color.yellow);

                bones[i].transform.localRotation = Quaternion.Euler(startAngle + angle1, 0, 0);
                //bones[i].transform.forward =  Quaternion.Euler(angle1, 0, 0) * dirToEnd;

                // law of sines
                float remainingLengthSquared = (end - bones[i].GetTail()).sqrMagnitude;
                float angle2 = 180 - Mathf.Rad2Deg * Mathf.Acos((bones[i].length * bones[i].length - distToEnd * distToEnd +
                    remainingLengthSquared) / (2.0f * bones[i].length * Mathf.Sqrt(remainingLengthSquared)));
                AngleRegion region2 = new AngleRegion(bones[i + 1].xLimit);
                AngleRegion angleRegion2 = new AngleRegion(-angle2, angle2);
                region2.Intersect(angleRegion2);
                //angle2 = region2.ClampInside(angle2);
                DebugRegion(i, region2, Color.yellow);
                DebugRegion(i, angleRegion2, Color.yellow);

                //Debug.Log(angle1 + " " + angle2);
                bones[i + 1].transform.forward = Quaternion.Euler(angle1 - Mathf.Sign(angle1) * angle2, 0, 0) * dirToEnd;
            }
        }

        private Vector3 PickWithContraints(int i)
        {
            Vector3 point = Vector3.forward;
            point = Quaternion.Euler(0, Random.Range(bones[i].yLimit.min, bones[i].yLimit.max), 0) * point;
            point = Quaternion.Euler(Random.Range(bones[i].xLimit.min, bones[i].xLimit.max), 0, 0) * point;
            return point;
        }

        private float GetLengthFrom(int i)
        {
            float length = 0;
            for (int j = i + 1; j < bones.Length; j++)
            {
                length += bones[j].length;
            }
            return length;
        }
        private float GetDistToEnd(int i)
        {
            return Mathf.Sqrt(GetDistToEndSquared(i));
        }
        private float GetDistToEndSquared(int i)
        {
            return (bones[i].transform.position - endEffector.transform.position).sqrMagnitude;
        }

        private void DebugRegion(int i, AngleRegion region, Color color)
        {
            if (i == -1)
            {
                region.Debug(bones[i + 1].transform.position, bones[i + 1].transform.parent.forward, bones[i + 1].transform.parent.right, color);
            }
            else
            {
                region.Debug(bones[i + 1].transform.position, bones[i].transform.forward, bones[i].transform.right, color);
            }
        }

        /*private void Resolve(int i, float prevAngle)
        {
            if (i < bones.Length - 2)
            {
                // law of cosines find alpha 
                float endDistSquared = GetEndDistSquared(i);
                float endDist = Mathf.Sqrt(endDistSquared);
                float remainingLength = GetRemainingLength(i);
                Debug.Log(i + ": " + endDist + " " + remainingLength);
                if (endDist > remainingLength + bones[i].length)
                {
                    bones[i].transform.forward = endEffector.transform.position - bones[i].transform.position;
                    ResolveExtended(i + 1);
                }
                else
                {
                    float angle = Mathf.Rad2Deg *
                        Mathf.Acos((GetLengthSquared(i) + endDistSquared - Mathf.Pow(remainingLength, 2)) /
                        (2.0f * bones[i].length * endDist));

                    if (i == 0)
                    {
                        angle = -angle * weight;
                        bones[i].transform.localRotation = Quaternion.Euler(angle, 0, 0);
                        Debug.DrawLine(bones[i].transform.position + bones[i].transform.forward * bones[i].length, endEffector.transform.position);
                        remainingLength = Vector3.Distance(bones[i].transform.position + bones[i].transform.forward * bones[i].length, endEffector.transform.position);
                        prevAngle = Mathf.Rad2Deg * Mathf.Asin(endDist * Mathf.Sin(Mathf.Deg2Rad * angle) / remainingLength);
                        Debug.Log(angle + " " + prevAngle);
                        Resolve(i + 1, -prevAngle);
                    }
                    else
                    {
                        angle =  prevAngle + angle * weight;
                        //bones[i].transform.forward = endEffector.transform.position - bones[i].transform.position;
                        //bones[i].transform.localRotation *= Quaternion.Euler(angle, 0, 0);
                        bones[i].transform.localRotation = Quaternion.Euler(angle, 0, 0);
                        Resolve(i + 1, angle);
                    }
                }
            }
            else
            {
                Resolve2Chain(i + 1);
            }
        }
        private void ResolveExtended(int i)
        {
            if (i < bones.Length)
            {
                bones[i].transform.localRotation = Quaternion.identity;
                ResolveExtended(i + 1);
            }
        }
        private void Resolve2Chain(int i)
        {

        }

        private float GetLengthSquared(int i)
        {
            return bones[i].length * bones[i].length;
        }
        private float GetRemainingLength(int i)
        {
            float length = 0;
            for (int j = i + 1; j < bones.Length; j++)
            {
                length += bones[j].length;
            }
            return length;
        }
        private float GetRemainingLengthSquared(int i) 
        {
            return Mathf.Pow(GetRemainingLength(i), 2);
        }
        private float GetEndDist(int i)
        {
            return Mathf.Sqrt(GetEndDistSquared(i));
            //return Mathf.Clamp((bones[i].transform.position - endEffector.transform.position).magnitude, ERROR, GetRemainingLength(i) - ERROR);
        }
        private float GetEndDistSquared(int i)
        {
            return (bones[i].transform.position - endEffector.transform.position).sqrMagnitude;
        }*/
    }
}
