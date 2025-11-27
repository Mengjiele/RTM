#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int isrc;
    float xsrc, ysrc, zsrc;
    float xrec, yrec, zrec;
	
	int iy_move, ishot;
	
    int iys, ixs;
    int iyr, ixr;
    
    // 配置参数
    int yshots = 8;       //观测系统的炮数
	int y_moves =15;
	int nys = yshots*y_moves;	//y方向的炮数
    int nxs = 60;     // 炮线方向需要移动115次
	
    int nxr = 320;      // 每条线2000个检波器
    int nyr = 40;    // 40条检波器线
    
    float dxr = 25;   // 检波器点间距(m)
    float dyr = 200;     // 检波器线间距(m)
    float oxr = 0;   // 检波器起始x位置(最小x坐标)
    float oyr = 0;   // 检波器起始y位置(最小y坐标)
    
    float dxs = 200;    // 炮线间距(m)
    float dys = 25;   // 炮点间距(m)
    float oys = oyr + (nyr/2)*dyr - (yshots/2)*dys; // 炮点起始y位置(居中)
    float oxs = oxr + (nxr/2)*dxr;   // 炮点起始x位置
	
	float y_move_step = 200;  // y方向移动步长(m)
	
    
    zsrc = 10;          // 震源深度位置(m)
    zrec = 10;         // 检波器深度位置(m)
    
    int ns = nxs * nys; // 总炮数
    int nr = nxr * nyr; // 每次放炮对应的检波器数
    
    printf("total shots (ns) = %d\n", ns);
    printf("geoghones numbers (ng) = %d\n", nr);
    
    // 打开震源文件
    FILE *src_file = fopen("sources.txt", "w");
    if (!src_file)
	{
        printf("error: source file cannot be created\n");
        return 1;
    }
    
    // 打开检波器文件
    FILE *rec_file = fopen("receivers.txt", "w");
    if (!rec_file)
	{
        printf("error: receivers file cannot be created\n");
        fclose(src_file);
        return 1;
    }
	
	fprintf(src_file, "shot_num\tz/m\t\tx/m\t\ty/m\n");
	fprintf(rec_file, "shot_num\tz/m\t\tx/m\t\ty/m\n");
	
	int shot_id = 0;
	for(ixs = 0; ixs < nxs; ixs++)
	{
		for (iy_move = 0; iy_move < y_moves; iy_move++)
		{
			float cur_y_offset = iy_move * y_move_step;
			
			float cur_xs_offset = oxs + ixs * dxs;
			float cur_xr_offset = ixs * dxs;
			
			// 计算当前炮点起始位置（位于阵列中央）
			float current_shot_start =	oys + cur_y_offset;
			
			// 计算当前检波器阵列起始位置
            float cur_oyr = oyr + cur_y_offset;
            float cur_oxr = oxr + cur_xr_offset;
			
			// 放当前排列的8炮
            for (ishot = 0; ishot < yshots; ishot++)
			{
				 // 炮点坐标
                ysrc = current_shot_start + ishot * dys;
                xsrc = cur_xs_offset;
				
                fprintf(src_file, "%-5d %5.2f %15.2f %15.2f\n", shot_id, zsrc, xsrc, ysrc);// 写入炮点信息
				
                for (int ixr = 0; ixr < nxr; ixr++)
				{
                    for (int iyr = 0; iyr < nyr; iyr++)
					{
                        xrec = cur_oxr + ixr * dxr;
                        yrec = cur_oyr + iyr * dyr;
                        
                        fprintf(rec_file, "%d %.2f %.2f %.2f\n", shot_id, zrec, xrec, yrec);// 写入对应的检波器信息
                    }
                }
                shot_id++;
			}
		}
	}
    
    
 /*   for (iys = 1; iys <= nys; iys++) // 遍历所有炮点位置(炮线方向)
	{
        ysrc = oys + (iys - 1) * dys;
		
        for (ixs = 1; ixs <= nxs; ixs++)// 遍历每条炮线上的8炮
		{
            xsrc = oxs + (ixs - 1) * dxs;
            isrc = ixs + nxs * (iys - 1);
            
            // 写入炮点信息
            fprintf(src_file, "%d %.2f %.2f %.2f\n", isrc, zsrc, xsrc, ysrc);
            
            for (ixr = 1; ixr <= nxr; ixr++)// 遍历所有检波器位置(固定排列)
			{
                xrec = oxr + (ixr - 1) * dxr;
                
                for (iyr = 1; iyr <= nyr; iyr++)
				{
                    yrec = oyr + (iyr - 1) * dyr;
                    
                    // 写入检波器信息(与当前炮对应)
                    fprintf(rec_file, "%d %.2f %.2f %.2f\n", isrc, zrec, xrec, yrec);
                }
            }
			
        }
    }*/
    
    fclose(src_file);
    fclose(rec_file);
    
    // 计算最大炮检距
    float max_x = oxr + (nxr-1)*dxr;
    float max_y = oyr + (nyr-1)*dyr;
    float max_offset = sqrt(pow(max_x - oxs, 2) + pow(max_y - oys, 2));
    
    printf("acquisition done:\n");
    printf("  - sources.txt: %d shot\n", ns);
    printf("  - receivers.txt: %d geophones\n", nr);
    printf("  - max_offset: %.2f m\n", max_offset);
    printf("  - coordinate range: x[%.2f, %.2f], y[%.2f, %.2f]\n", 
           fmin(oxs, oxr), fmax(oxs+(nxs-1)*dxs, oxr+(nxr-1)*dxr),
           fmin(oys, oyr), fmax(oys+(nys-1)*dys, oyr+(nyr-1)*dyr));
    
    return 0;
}
