using UnityEngine;
using InverseKinematics2D;

public class LizardLeg : MonoBehaviour
{
    public Transform target;
    public bool flipSide;
    
    public IK2Segment start;
    public IK2Segment end;
    public IK2Segment foot;

    private Vector2 prevPos = Vector2.zero;
    private float prevRot = 0;

    private Lizard lizard;

    private Vector2 deltaP;
    private float deltaR;

    private void Awake()
    {
        lizard = GetComponentInParent<Lizard>();
        prevPos = transform.position;
        prevRot = transform.eulerAngles.z;
        target.position = RestPos;
    }

    public bool debv;
    private void Update()
    {
        float dist = Mathf.Clamp(Vector3.Distance(transform.position, target.position), float.Epsilon, MaxDist);
        float angleToTarget = Vector3.SignedAngle(target.position - transform.position, transform.right, transform.forward);

        float a1 = Mathf.Acos((end.length * end.length + dist * dist - start.length * start.length) / 
            (2 * end.length * dist)) * Mathf.Rad2Deg;
        float a2 = 180 + Mathf.Acos((end.length * end.length + start.length * start.length - dist * dist) /
            (2 * end.length * start.length)) * Mathf.Rad2Deg;

        start.transform.localRotation = Quaternion.Euler(0, 0, !flipSide ? -a1 - angleToTarget : a1 - angleToTarget);
        end.transform.localRotation = Quaternion.Euler(0, 0, !flipSide ? -a2 : +a2);

        deltaP = (Vector2)transform.position - prevPos;
        deltaR = transform.eulerAngles.z - prevRot;

        debv = OutOfRange;
        prevPos = transform.position;
        prevRot = transform.eulerAngles.z;
    }

    public float MaxDist => start.length + end.length;
    public Vector2 DeltaPos => deltaP;
    public float DeltaRot => deltaR;
    public Vector2 RestPos
    {
        get
        {
            Vector2 v1 = IK2.Vector(start.IdealAngle, start.length);
            Vector2 v2 = IK2.Vector(start.IdealAngle + end.IdealAngle, end.length);
            Vector2 local = Quaternion.Euler(0, 0, DeltaRot) * (v1 + v2);
            return (Vector2)transform.TransformPoint(local) + DeltaPos.normalized / 2;
        }
    }
    public bool OutOfRange => Vector2.Distance(RestPos, target.position) > 1.5f;

    private void OnDrawGizmos()
    {
        if (lizard != null)
        {
            Gizmos.DrawSphere(RestPos, 0.5f);
        }
    }
}
