using UnityEngine;
using InverseKinematicsPolarTest;

public class CircleIntersectionTest : MonoBehaviour
{
    public Transform circle1;
    public Transform circle2;
    public float r1 = 1;
    public float r2 = 1;

    void Update()
    {
        Debug.DrawCircle(circle1.position, r1, 128, Color.red);
        Debug.DrawCircle(circle2.position, r2, 128, Color.green);

        CircleIntersection intersection = IKUtility.TwoCircleIntersection(circle1.position, circle2.position, r1, r2);
        Debug.DrawPoints(intersection, Color.blue);
    }
}
