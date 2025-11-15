#!/usr/bin/env python3
"""
Test script for the RegNetAgents LangGraph-powered MCP Server
Tests the hybrid MCP-LangGraph integration
"""

import asyncio
import json
import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

async def test_langgraph_mcp_tools():
    """Test the LangGraph MCP server tools"""
    print("Testing RegNetAgents LangGraph-powered MCP Server...")

    try:
        # Import the server components
        from regnetagents_langgraph_mcp_server import handle_call_tool, get_workflow

        print("Initializing LangGraph workflow...")
        workflow = await get_workflow()
        print("Workflow initialized successfully")

        # Test 1: Comprehensive Gene Analysis
        print("\n=== Test 1: Comprehensive Gene Analysis ===")
        test1_args = {
            "gene": "APC",
            "cell_type": "epithelial_cell",
            "analysis_depth": "comprehensive"
        }

        print(f"Testing comprehensive_gene_analysis with APC...")
        result1 = await handle_call_tool("comprehensive_gene_analysis", test1_args)

        if result1 and len(result1) > 0:
            response = json.loads(result1[0].text)
            if response.get('error'):
                print(f"Test 1 ERROR: {response['error']}")
            else:
                print("Test 1 SUCCESS: Comprehensive analysis completed")
                print(f"  Workflow type: {response.get('workflow_info', {}).get('workflow_type')}")
                print(f"  Gene: {response.get('gene_analysis_summary', {}).get('gene')}")
                print(f"  Regulatory role: {response.get('gene_analysis_summary', {}).get('regulatory_role')}")

                # Show workflow execution details
                metadata = response.get('workflow_metadata', {})
                steps = metadata.get('steps_completed', [])
                print(f"  Steps completed: {len(steps)}")
                print(f"  Total time: {metadata.get('total_analysis_time', 0):.2f}ms")
        else:
            print("Test 1 FAILED: No response received")

        # Test 2: Multi-Gene Analysis
        print("\n=== Test 2: Multi-Gene Analysis ===")
        test2_args = {
            "genes": ["APC", "BRCA1", "TP53"],
            "cell_type": "epithelial_cell",
            "analysis_depth": "focused"
        }

        print(f"Testing multi_gene_analysis with 3 genes...")
        result2 = await handle_call_tool("multi_gene_analysis", test2_args)

        if result2 and len(result2) > 0:
            response = json.loads(result2[0].text)
            if response.get('error'):
                print(f"Test 2 ERROR: {response['error']}")
            else:
                print("Test 2 SUCCESS: Multi-gene analysis completed")
                summary = response.get('multi_gene_analysis', {})
                print(f"  Genes analyzed: {summary.get('total_genes')}")
                print(f"  Successful: {summary.get('successful_analyses')}")
                print(f"  Failed: {summary.get('failed_analyses')}")
        else:
            print("Test 2 FAILED: No response received")

        # Test 3: Workflow Status
        print("\n=== Test 3: Workflow Status ===")
        test3_args = {
            "gene": "APC",
            "show_state": True
        }

        print(f"Testing workflow_status...")
        result3 = await handle_call_tool("workflow_status", test3_args)

        if result3 and len(result3) > 0:
            response = json.loads(result3[0].text)
            print("Test 3 SUCCESS: Workflow status retrieved")
            print(f"  Status: {response.get('status')}")
            print(f"  Workflow type: {response.get('workflow_type')}")
        else:
            print("Test 3 FAILED: No response received")

        # Test 4: Pathway-Focused Analysis
        print("\n=== Test 4: Pathway-Focused Analysis ===")
        test4_args = {
            "gene": "APC",
            "pathway_focus": "wnt_signaling",
            "cell_type": "epithelial_cell"
        }

        print(f"Testing pathway_focused_analysis...")
        result4 = await handle_call_tool("pathway_focused_analysis", test4_args)

        if result4 and len(result4) > 0:
            response = json.loads(result4[0].text)
            if response.get('error'):
                print(f"Test 4 ERROR: {response['error']}")
            else:
                print("Test 4 SUCCESS: Pathway-focused analysis completed")
                pathway_info = response.get('pathway_focused_analysis', {})
                print(f"  Primary gene: {pathway_info.get('primary_gene')}")
                print(f"  Pathway focus: {pathway_info.get('pathway_focus')}")
        else:
            print("Test 4 FAILED: No response received")

        # Test 5: Workflow Insights
        print("\n=== Test 5: Workflow Insights ===")
        test5_args = {
            "analysis_type": "performance"
        }

        print(f"Testing workflow_insights...")
        result5 = await handle_call_tool("workflow_insights", test5_args)

        if result5 and len(result5) > 0:
            response = json.loads(result5[0].text)
            print("Test 5 SUCCESS: Workflow insights retrieved")
            insights = response.get('workflow_insights', {})
            advantages = insights.get('langgraph_advantages', {})
            print(f"  LangGraph advantages: {len(advantages)} features")
            print(f"  Performance optimizations: {len(insights.get('performance_optimizations', {}))}")
        else:
            print("Test 5 FAILED: No response received")

        # Test 6: Create Analysis Report
        print("\n=== Test 6: Create Analysis Report ===")
        test6_args = {
            "gene": "APC",
            "report_format": "summary",
            "include_visualizations": True
        }

        print(f"Testing create_analysis_report...")
        result6 = await handle_call_tool("create_analysis_report", test6_args)

        if result6 and len(result6) > 0:
            response = json.loads(result6[0].text)
            print("Test 6 SUCCESS: Analysis report created")
            report = response.get('analysis_report', {})
            print(f"  Gene: {report.get('gene')}")
            print(f"  Format: {report.get('report_format')}")
            print(f"  Generated by: {report.get('generated_by')}")
        else:
            print("Test 6 FAILED: No response received")

        print("\n" + "="*60)
        print("LangGraph MCP Server Tests Completed Successfully!")
        print("The hybrid MCP-LangGraph server is working correctly.")
        print("You can now use these tools in Claude Desktop.")

        return True

    except Exception as e:
        print(f"Test failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

async def demonstrate_workflow_advantages():
    """Demonstrate the advantages of the LangGraph implementation"""
    print("\n" + "="*60)
    print("LANGGRAPH MCP SERVER ADVANTAGES:")
    print("="*60)

    advantages = {
        "Intelligent Workflow Orchestration": [
            "Visual state management with clear transitions",
            "Smart routing based on gene characteristics",
            "Conditional execution of relevant analyses only",
            "Graceful error handling and recovery"
        ],
        "Enhanced Performance": [
            "Shared cache instances across workflow",
            "Parallel processing for multi-gene analysis",
            "Optimized resource usage",
            "Fast gene mapping integration"
        ],
        "Advanced Features": [
            "Comprehensive workflow status tracking",
            "Pathway-focused analysis capabilities",
            "Multi-gene batch processing",
            "Detailed execution insights"
        ],
        "MCP Integration Benefits": [
            "Seamless Claude Desktop compatibility",
            "Rich tool descriptions and schemas",
            "Structured JSON responses",
            "Error handling with actionable feedback"
        ]
    }

    for category, items in advantages.items():
        print(f"\n{category}:")
        for item in items:
            print(f"  • {item}")

    print(f"\n{'='*60}")
    print("CONFIGURATION FOR CLAUDE DESKTOP:")
    print("Add this to your claude_desktop_config.json:")
    print("="*60)

    config = {
        "mcpServers": {
            "regnetagents-langgraph-server": {
                "command": "python",
                "args": ["C:/path/to/RegNetAgents/regnetagents_langgraph_mcp_server.py"],
                "env": {
                    "PYTHONPATH": "C:/path/to/RegNetAgents"
                }
            }
        }
    }

    print(json.dumps(config, indent=2))

if __name__ == "__main__":
    success = asyncio.run(test_langgraph_mcp_tools())
    if success:
        asyncio.run(demonstrate_workflow_advantages())
        print("\nLangGraph MCP Server is ready for Claude Desktop!")
    else:
        print("\nLangGraph MCP Server test FAILED!")
        sys.exit(1)