using UnityEngine;

public class Lizard : MonoBehaviour
{
    public enum LegMode { Ground, Air }
    public LizardLeg[] legs;
    private LegMode[] modes;

    public int groundCount = 2;

    public float legMoveSpeed = 3f;
    public float moveSpeed = 1f;
    public float rotateSpeed = 1f;

    private void Awake()
    {
        modes = new LegMode[legs.Length];
    }

    private void Update()
    {
        Vector2 mousePos = Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z));
        Vector3 faceDir = mousePos - (Vector2)transform.position;
        // move and rotate lizard toward target point, mouse
        if (Input.GetKey(KeyCode.Space))
        {
            transform.Translate(Vector2.right * moveSpeed * Time.deltaTime, Space.Self);
            float faceDirAngle = Vector3.SignedAngle(transform.right, faceDir, Vector3.back);
            faceDirAngle = Mathf.Clamp(faceDirAngle, -rotateSpeed, rotateSpeed);
            transform.Rotate(Vector3.back, faceDirAngle * Time.deltaTime, Space.Self);
        }

        for (int i = 0; i < legs.Length; i++)
        {
            if (modes[i] == LegMode.Ground)
            {
                if (legs[i].OutOfRange && LegsOnGround > groundCount)
                {
                    modes[i] = LegMode.Air;
                }
            }
            else if (modes[i] == LegMode.Air)
            {
                legs[i].target.position = Vector2.MoveTowards(legs[i].target.position, legs[i].RestPos, legMoveSpeed * Time.deltaTime);
                if (Vector3.Distance(legs[i].target.position, legs[i].RestPos) < 0.1f)
                {
                    modes[i] = LegMode.Ground;
                }
            }
        }
    }

    public int LegsOnGround
    {
        get
        {
            int count = 0;
            for (int i = 0; i < legs.Length; i++)
            {
                if (modes[i] == LegMode.Ground)
                {
                    count++;
                }
            }
            return count;
        }
    }
}
